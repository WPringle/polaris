from mpas_tools.io import open_dataset, write_netcdf

from polaris import Step


class CoarsenBathymetryStep(Step):
    """
    A step for block-averaging a lake's native-resolution bathymetry down
    to a coarser grid before the conservative remap.

    The native DEM (e.g. NOAA's 3 arc-second, ~90 m, Great Lakes grids) is
    far finer than this task's own finest target cell width
    (``coastal_resolution_km``, 925 m). Conservative-remap accuracy only
    needs the source to be modestly finer than the target, not ~10x
    finer, and the remap's cost (via ESMF/MOAB's ``mbtempest``) scales
    with source point count, which for a lake-sized DEM at native
    resolution can run into the tens of millions of points and take a
    very long time. Coarsening to ``coarsen_resolution_km`` (still several
    times finer than the coastal target) cuts that cost substantially
    while remaining fine enough for an area-weighted average.

    This does not affect the mesh step's own cell-width shaping or the
    cull step's wet-mask classification, both of which use nearest-
    neighbor/KD-tree lookups against the native-resolution bathymetry
    directly and are not the slow part.
    """

    def __init__(self, component, subdir, config, bathymetry_step, lake_name):
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
            The step that prepares the lake's native-resolution bathymetry

        lake_name : str
            The name of the lake, used only for log messages
        """
        super().__init__(component=component, name='coarsen', subdir=subdir)
        self.lake_name = lake_name
        self.set_shared_config(config, link='great_lakes_mesh.cfg')

        self.add_input_file(
            filename='bathymetry.nc',
            work_dir_target=f'{bathymetry_step.path}/bathymetry.nc',
        )
        self.add_output_file('bathymetry_coarse.nc')

    def run(self):
        """
        Block-average the bathymetry down to ``coarsen_factor``-times
        coarser resolution
        """
        logger = self.logger
        section = self.config['great_lakes_topo']
        factor = section.getint('coarsen_factor')

        if factor <= 1:
            logger.info(
                f'coarsen_factor <= 1; using native-resolution bathymetry '
                f'for {self.lake_name} unchanged'
            )
            with open_dataset('bathymetry.nc') as ds:
                write_netcdf(ds, 'bathymetry_coarse.nc')
            return

        logger.info(
            f'Coarsening {self.lake_name} bathymetry by a factor of {factor}'
        )
        with open_dataset('bathymetry.nc') as ds:
            ds_coarse = ds.coarsen(
                lon=factor, lat=factor, boundary='trim'
            ).mean()
            # coarsen().mean() does not preserve attributes
            ds_coarse.lon.attrs['units'] = 'degrees_east'
            ds_coarse.lat.attrs['units'] = 'degrees_north'
            write_netcdf(ds_coarse, 'bathymetry_coarse.nc')
