import threading
import time

from mpas_tools.io import open_dataset

from polaris.remap import MappingFileStep


class GreatLakesTopoMap(MappingFileStep):
    """
    A step for making a mapping file from a lake's lon/lat bathymetry grid
    to an MPAS mesh

    Attributes
    ----------
    mesh_name : str
        The name of the mesh
    """

    def __init__(
        self,
        component,
        name,
        subdir,
        config,
        mesh_name,
        bathymetry_step,
        mesh_step,
        mesh_filename,
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

        mesh_name : str
            The name of the mesh

        bathymetry_step : polaris.Step
            The step that provides the source bathymetry file

        mesh_step : polaris.Step
            The mesh step (base or culled) to map onto

        mesh_filename : str
            The mesh file name within ``mesh_step``'s work directory

        bathymetry_filename : str, optional
            The bathymetry file name within ``bathymetry_step``'s work
            directory. Pass the output of
            :py:class:`polaris.tasks.ocean.great_lakes.topo.coarsen.CoarsenBathymetryStep`
            (``bathymetry_coarse.nc``) rather than the native-resolution
            file for any remap, since the source grid size dominates
            ``mbtempest``'s runtime.
        """  # noqa: E501
        super().__init__(
            component=component,
            name=name,
            subdir=subdir,
            ntasks=1,
            min_tasks=1,
            method='conserve',
        )
        self.mesh_name = mesh_name
        self.set_shared_config(config, link='great_lakes_mesh.cfg')

        self.add_input_file(
            filename='bathymetry.nc',
            work_dir_target=f'{bathymetry_step.path}/{bathymetry_filename}',
        )
        self.add_input_file(
            filename='mesh.nc',
            work_dir_target=f'{mesh_step.path}/{mesh_filename}',
        )

    def run(self):
        """
        Set up the source and destination grids for this step
        """
        self.remapper.src_from_lon_lat(
            filename='bathymetry.nc',
            mesh_name='great_lakes_bathymetry',
            lon_var='lon',
            lat_var='lat',
        )
        self.remapper.dst_from_mpas(
            filename='mesh.nc', mesh_name=self.mesh_name
        )

        with open_dataset('bathymetry.nc') as ds_src:
            src_points = ds_src.sizes['lon'] * ds_src.sizes['lat']
        with open_dataset('mesh.nc') as ds_dst:
            dst_cells = ds_dst.sizes['nCells']
        self.logger.info(
            f'Computing conservative remap weights: {src_points} source '
            f'bathymetry points -> {dst_cells} destination mesh cells. '
            f'mbtempest reports no completion fraction, so progress below '
            f'is only an elapsed-time liveness check, not a percentage.'
        )

        interval_s = self.config.getfloat(
            'great_lakes_topo', 'remap_progress_interval_s'
        )
        stop = threading.Event()

        def _heartbeat():
            start = time.monotonic()
            while not stop.wait(interval_s):
                elapsed = time.monotonic() - start
                self.logger.info(
                    f'  ... still computing remap weights '
                    f'({elapsed:.0f} s elapsed; mbtempest gives no '
                    f'completion fraction, this is only a liveness check)'
                )

        heartbeat = threading.Thread(target=_heartbeat, daemon=True)
        heartbeat.start()
        try:
            super().run()
        finally:
            stop.set()
            heartbeat.join()
