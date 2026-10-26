import xarray as xr
from mpas_tools.io import open_dataset, write_netcdf

from polaris import Step


class GreatLakesTopoRemap(Step):
    """
    A step for remapping a lake's bathymetry onto an MPAS mesh

    Attributes
    ----------
    bathymetry_step : polaris.Step
        The step that prepares the lake's bathymetry

    topo_map : polaris.tasks.ocean.great_lakes.topo.map.GreatLakesTopoMap
        The step that creates the mapping file used to remap the
        bathymetry
    """

    def __init__(
        self,
        component,
        name,
        subdir,
        config,
        bathymetry_step,
        topo_map,
        bathymetry_filename='bathymetry.nc',
    ):
        """
        Create the step

        Parameters
        ----------
        component : polaris.Component
            The component the step belongs to

        name : str
            The name of the step

        subdir : str
            The subdirectory for the step

        config : polaris.config.PolarisConfigParser
            A shared config parser

        bathymetry_step : polaris.Step
            The step that prepares the lake's bathymetry

        topo_map : polaris.tasks.ocean.great_lakes.topo.map.GreatLakesTopoMap
            The step that creates the mapping file used to remap the
            bathymetry

        bathymetry_filename : str, optional
            The bathymetry file name within ``bathymetry_step``'s work
            directory. Should match whatever file ``topo_map`` was built
            from (e.g. the output of
            :py:class:`polaris.tasks.ocean.great_lakes.topo.coarsen.CoarsenBathymetryStep`)
            so the mapping file's source grid matches the data being
            remapped.
        """  # noqa: E501
        super().__init__(component=component, name=name, subdir=subdir)
        self.set_shared_config(config, link='great_lakes_mesh.cfg')
        self.bathymetry_step = bathymetry_step
        self.topo_map = topo_map

        self.add_input_file(
            filename='bathymetry.nc',
            work_dir_target=f'{bathymetry_step.path}/{bathymetry_filename}',
        )
        self.add_dependency(topo_map, name='topo_map')
        self.add_output_file('topography_remapped.nc')

    def run(self):
        """
        Run this step
        """
        self._preprocess()
        self._remap()
        self._rename_and_renormalize()

    @staticmethod
    def _preprocess():
        with xr.open_dataset('bathymetry.nc') as ds_in:
            ds_in = ds_in.rename({'elevation': 'base_elevation'})
            ocean_mask = (ds_in.base_elevation < 0.0).astype(float)
            ds_in['oceanFraction'] = ocean_mask
            # Zero out land's contribution before the conservative remap so
            # the remapped elevation only integrates over the wet area.
            # Land elevation near Great Lakes shorelines is tens of meters
            # positive (unlike isomip_plus's near-zero ice-shelf land), so
            # leaving it in would bias the post-remap division below toward
            # positive values for any partially-land cell.
            ds_in['base_elevation'] = ds_in.base_elevation * ocean_mask
            ds_in.to_netcdf('bathymetry_processed.nc')

    def _remap(self):
        logger = self.logger
        topo_map = self.dependencies['topo_map']
        remapper = topo_map.remapper

        remapper.ncremap(
            in_filename='bathymetry_processed.nc',
            out_filename='topography_ncremap.nc',
            overwrite=True,
            logger=logger,
        )

    @staticmethod
    def _rename_and_renormalize():
        with open_dataset('topography_ncremap.nc') as ds_out:
            drop_candidates = [
                'lon',
                'lat',
                'area',
                'lat_vertices',
                'lon_vertices',
            ]
            drop = [var for var in drop_candidates if var in ds_out]
            ds_out = ds_out.drop_vars(drop)

            ocean_frac = ds_out.oceanFraction
            mask = ocean_frac > 0.0
            attrs = ds_out.base_elevation.attrs
            ds_out['base_elevation'] = (
                ds_out.base_elevation / ocean_frac
            ).where(mask)
            ds_out.base_elevation.attrs = attrs

            write_netcdf(ds_out, 'topography_remapped.nc')
