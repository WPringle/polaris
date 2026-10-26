import numpy as np
import xarray as xr
from mpas_tools.io import open_dataset, write_netcdf
from mpas_tools.logging import LoggingContext
from mpas_tools.mesh.conversion import cull
from mpas_tools.mesh.creation.sort_mesh import sort_mesh

from polaris import Step
from polaris.model_step import make_graph_file


class GreatLakesCullMesh(Step):
    """
    A step for culling a lake's base mesh down to its wet cells

    Attributes
    ----------
    base_mesh : polaris.Step
        The base mesh step containing input files to this step

    wet_mask_step : polaris.Step
        The step that provides a point-sampled wet/dry mask on the base
        mesh, used for the cull decision
    """

    def __init__(self, component, subdir, config, base_mesh, wet_mask_step):
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

        base_mesh : polaris.Step
            The base mesh step containing input files to this step

        wet_mask_step : polaris.Step
            The step that provides a point-sampled wet/dry mask on the
            base mesh, used for the cull decision
        """
        super().__init__(component=component, name='cull_mesh', subdir=subdir)
        self.base_mesh = base_mesh
        self.wet_mask_step = wet_mask_step
        self.set_shared_config(config, link='great_lakes_mesh.cfg')

        self.add_input_file(
            filename='base_mesh.nc',
            work_dir_target=f'{base_mesh.path}/base_mesh.nc',
        )
        self.add_input_file(
            filename='topography.nc',
            work_dir_target=f'{wet_mask_step.path}/ocean_fraction.nc',
        )

        for file in ['culled_mesh.nc', 'culled_graph.info']:
            self.add_output_file(filename=file)

    def run(self):
        """
        Run this step
        """
        logger = self.logger
        section = self.config['great_lakes_topo']
        min_ocean_fraction = section.getfloat('min_ocean_fraction')

        with LoggingContext(name=__name__, logger=logger) as logger:
            _land_mask_from_topo(
                topo_filename='topography.nc',
                mask_filename='land_mask.nc',
                min_ocean_fraction=min_ocean_fraction,
            )

            ds_base = open_dataset('base_mesh.nc')
            ds_land_mask = open_dataset('land_mask.nc')

            # cull the mesh based on the land mask
            ds_culled = cull(
                ds_base, dsMask=ds_land_mask, logger=logger, dir='.'
            )

            # sort the cell, edge and vertex indices for better performance
            ds_culled = sort_mesh(ds_culled)

            out_filename = 'culled_mesh.nc'
            write_netcdf(ds_culled, out_filename)

            # we need to make the graph file after sorting
            make_graph_file(
                mesh_filename='culled_mesh.nc',
                graph_filename='culled_graph.info',
            )


def _land_mask_from_topo(topo_filename, mask_filename, min_ocean_fraction):
    ds_topo = open_dataset(topo_filename)

    ocean_frac = ds_topo.oceanFraction

    # we want the mask to be 1 where there's not ocean
    mask = np.logical_or(ocean_frac < min_ocean_fraction, ocean_frac.isnull())
    cull_mask = xr.where(mask, 1, 0)

    cull_mask = cull_mask.expand_dims(dim='nRegions', axis=1)

    ds_mask = xr.Dataset()
    ds_mask['regionCellMasks'] = cull_mask
    write_netcdf(ds_mask, mask_filename)
