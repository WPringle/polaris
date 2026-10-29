import gzip
import shutil

import numpy as np
import xarray as xr
from mpas_tools.io import write_netcdf

from polaris import Step


class PrepareBathymetryStep(Step):
    """
    A step for preparing a Great Lakes basin's bathymetry: decompressing
    the source file(s) if needed, standardizing them to ``lon``/``lat``/
    ``elevation`` names, and merging multiple sources onto one grid for a
    basin that spans more than one downloaded file (see
    :py:data:`polaris.tasks.ocean.great_lakes.lakes.GREAT_LAKES`).

    Attributes
    ----------
    lake_name : str
        The name of the lake, used only for log messages
    """

    def __init__(self, component, subdir, config, lake_name):
        """
        Create the step

        Parameters
        ----------
        component : polaris.Component
            The component the step belongs to

        subdir : str
            The subdirectory for the step

        config : polaris.config.PolarisConfigParser
            A shared config parser

        lake_name : str
            The name of the lake, used only for log messages
        """
        super().__init__(component=component, name='bathymetry', subdir=subdir)
        self.lake_name = lake_name
        self.set_shared_config(config, link='great_lakes_mesh.cfg')
        self.add_output_file('bathymetry.nc')

    def run(self):
        """
        Decompress (if needed), standardize, and (if more than one source
        is given) merge the basin's bathymetry file(s)

        The source path(s) come from the ``bathymetry_filenames`` config
        option, which is only known once the user has edited the work-dir
        config (it defaults to the ``<<<missing>>>`` sentinel), so it is
        read here rather than at ``setup()`` time -- the same pattern used
        by ``ComputeFeatureMasksStep`` for its ``mesh_filename`` option.
        """
        logger = self.logger
        section = self.config['great_lakes_topo']
        bathymetry_filenames = _require_filenames(
            section.getlist('bathymetry_filenames', fallback=[])
        )
        elevation_var = section.get('elevation_var')
        lon_var = section.get('lon_var')
        lat_var = section.get('lat_var')

        datasets = []
        for index, bathymetry_filename in enumerate(bathymetry_filenames):
            datasets.append(
                _read_one_source(
                    bathymetry_filename=bathymetry_filename,
                    elevation_var=elevation_var,
                    lon_var=lon_var,
                    lat_var=lat_var,
                    decompressed_filename=f'bathymetry_decompressed_{index}.grd',
                    lake_name=self.lake_name,
                    logger=logger,
                )
            )

        if len(datasets) > 1:
            logger.info(
                f'Merging {len(datasets)} bathymetry sources for '
                f'{self.lake_name}'
            )
            ds = _merge_bathymetry_grids(datasets)
        else:
            ds = datasets[0]

        lon_max_raw = (section.get('lon_max_deg', fallback='') or '').strip()
        if lon_max_raw:
            lon_max_deg = float(lon_max_raw)
            n_before = ds.sizes['lon']
            ds = ds.sel(lon=slice(None, lon_max_deg))
            logger.info(
                f'Clipping {self.lake_name} bathymetry to lon <= '
                f'{lon_max_deg} deg: dropped {n_before - ds.sizes["lon"]} '
                f'of {n_before} longitude columns'
            )

        # pyremap's LatLonGridDescriptor requires CF-compliant units on
        # the coordinate variables; GMT grids don't set them
        ds.lon.attrs['units'] = 'degrees_east'
        ds.lat.attrs['units'] = 'degrees_north'

        write_netcdf(ds, 'bathymetry.nc')


def _read_one_source(
    bathymetry_filename,
    elevation_var,
    lon_var,
    lat_var,
    decompressed_filename,
    lake_name,
    logger,
):
    with open(bathymetry_filename, 'rb') as f:
        is_gzipped = f.read(2) == b'\x1f\x8b'

    if is_gzipped:
        logger.info(
            f'Decompressing {lake_name} bathymetry: {bathymetry_filename}'
        )
        with gzip.open(bathymetry_filename, 'rb') as f_in:
            with open(decompressed_filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        decompressed_filename = bathymetry_filename

    with xr.open_dataset(decompressed_filename) as ds:
        rename = {}
        if elevation_var != 'elevation':
            rename[elevation_var] = 'elevation'
        if lon_var != 'lon':
            rename[lon_var] = 'lon'
        if lat_var != 'lat':
            rename[lat_var] = 'lat'
        if rename:
            ds = ds.rename(rename)

        ds = ds[['elevation']]
        ds = ds.transpose('lat', 'lon')
        ds.load()

    return ds


def _merge_bathymetry_grids(datasets):
    """
    Merge one or more lon/lat bathymetry grids onto a single grid.

    The datasets must all share the same lon/lat spacing (as NOAA's
    per-lake Great Lakes grids do), but may cover different, only
    partially overlapping extents -- e.g. Lake Michigan's and Lake
    Huron's grids, which tile together exactly at the Straits of
    Mackinac, down to an identical shared boundary column. Returns a
    single dataset spanning the union of all input extents, with NaN
    wherever none of the inputs has data (already handled downstream:
    cell classification checks ``np.isfinite`` before testing
    ``elevation < 0``, so NaN is treated the same as land).
    """
    dlon = float(datasets[0].lon[1] - datasets[0].lon[0])
    dlat = float(datasets[0].lat[1] - datasets[0].lat[0])
    for ds in datasets[1:]:
        if not np.isclose(
            float(ds.lon[1] - ds.lon[0]), dlon
        ) or not np.isclose(float(ds.lat[1] - ds.lat[0]), dlat):
            raise ValueError(
                'All bathymetry sources for a combined basin must share '
                'the same grid spacing.'
            )

    lon_min = min(float(ds.lon.min()) for ds in datasets)
    lon_max = max(float(ds.lon.max()) for ds in datasets)
    lat_min = min(float(ds.lat.min()) for ds in datasets)
    lat_max = max(float(ds.lat.max()) for ds in datasets)

    nlon = int(round((lon_max - lon_min) / dlon)) + 1
    nlat = int(round((lat_max - lat_min) / dlat)) + 1
    lon = lon_min + dlon * np.arange(nlon)
    lat = lat_min + dlat * np.arange(nlat)

    elevation = np.full((nlat, nlon), np.nan, dtype=np.float32)
    for ds in datasets:
        lon_start = int(round((float(ds.lon.min()) - lon_min) / dlon))
        lat_start = int(round((float(ds.lat.min()) - lat_min) / dlat))
        n_lat_ds, n_lon_ds = ds.elevation.shape
        elevation[
            lat_start : lat_start + n_lat_ds,
            lon_start : lon_start + n_lon_ds,
        ] = ds.elevation.values

    return xr.Dataset(
        {'elevation': (('lat', 'lon'), elevation)},
        coords={'lat': lat, 'lon': lon},
    )


def _require_filenames(filenames):
    if not filenames or filenames == ['<<<missing>>>']:
        raise ValueError(
            '[great_lakes_topo] bathymetry_filenames must be set to the '
            'path of one or more downloaded Great Lakes bathymetry files '
            'before this step can run.'
        )
    return filenames
