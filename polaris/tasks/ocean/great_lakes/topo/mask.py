import numpy as np
import xarray as xr
from mpas_tools.io import open_dataset, write_netcdf

from polaris import Step
from polaris.mesh.connectivity import connected_to_seeds, seed_mask_from_points


class WetMaskOnMeshStep(Step):
    """
    A step for building a cheap, point-sampled wet/dry mask on the
    uncull base mesh, used only to decide which cells
    :py:class:`polaris.tasks.ocean.great_lakes.cull.GreatLakesCullMesh`
    keeps.

    A full conservative remap (as used for the final ``base_elevation``/
    ``oceanFraction`` fields on the culled mesh) is not needed here: the
    cull decision only needs to know roughly which cells are wet, not an
    area-accurate ocean fraction, and the base mesh has far more cells
    than will survive culling. Point-sampling the bathymetry at each
    cell's center avoids running the expensive ESMF/MOAB conservative
    remap against the bathymetry's full native grid (potentially tens of
    millions of points) on a mesh that is mostly discarded anyway.
    """

    def __init__(self, component, subdir, config, bathymetry_step, base_mesh):
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

        bathymetry_step : polaris.Step
            The step that prepares the lake's bathymetry

        base_mesh : polaris.Step
            The uncull base mesh step
        """
        super().__init__(component=component, name='wet_mask', subdir=subdir)
        self.set_shared_config(config, link='great_lakes_mesh.cfg')

        self.add_input_file(
            filename='bathymetry.nc',
            work_dir_target=f'{bathymetry_step.path}/bathymetry.nc',
        )
        self.add_input_file(
            filename='base_mesh.nc',
            work_dir_target=f'{base_mesh.path}/base_mesh.nc',
        )
        self.add_output_file('ocean_fraction.nc')

    def run(self):
        """
        Run this step
        """
        with (
            open_dataset('bathymetry.nc') as ds_bathy,
            open_dataset('base_mesh.nc') as ds_mesh,
        ):
            lon_cell = np.degrees(ds_mesh.lonCell.values)
            # MPAS lonCell is in [0, 2*pi); shift to (-180, 180] to match
            # the bathymetry's longitude convention
            lon_cell = np.where(lon_cell > 180.0, lon_cell - 360.0, lon_cell)
            lat_cell = np.degrees(ds_mesh.latCell.values)

            elevation = ds_bathy.elevation.values
            bathy_lon = ds_bathy.lon.values
            bathy_lat = ds_bathy.lat.values

            sampled = ds_bathy.elevation.sel(
                lon=xr.DataArray(lon_cell, dims='nCells'),
                lat=xr.DataArray(lat_cell, dims='nCells'),
                method='nearest',
            ).values

            # ``.sel(..., method='nearest')`` always returns the closest
            # available point, even for a cell far outside the
            # bathymetry's own extent (e.g. anywhere in the coarse global
            # background mesh, which is most of the base mesh). Without
            # this bound, a distant cell could be nearest-matched to a wet
            # pixel at the edge of the lake's own small grid and
            # spuriously marked ocean.
            in_extent = (
                (lon_cell >= bathy_lon.min())
                & (lon_cell <= bathy_lon.max())
                & (lat_cell >= bathy_lat.min())
                & (lat_cell <= bathy_lat.max())
            )
            cell_mask = in_extent & np.isfinite(sampled) & (sampled < 0.0)

            # a coarse cell near the graded grid's transition zone can
            # have its single center point land on water while none of
            # its neighbors do, leaving an isolated "wet" cell that
            # breaks cull()/sort_mesh() downstream. Keep only cells
            # connected to the lake's deepest point, which is
            # unambiguously wet.
            deepest = np.unravel_index(
                np.nanargmin(elevation), elevation.shape
            )
            lon_seed = bathy_lon[deepest[1]]
            lat_seed = bathy_lat[deepest[0]]
            seed_mask = seed_mask_from_points(ds_mesh, lon_seed, lat_seed)
            connected = connected_to_seeds(ds_mesh, cell_mask, seed_mask)

        ocean_fraction = np.where(connected, 1.0, 0.0)

        ds_out = xr.Dataset()
        ds_out['oceanFraction'] = ('nCells', ocean_fraction)
        write_netcdf(ds_out, 'ocean_fraction.nc')
