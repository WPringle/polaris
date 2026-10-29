import numpy as np
import xarray as xr
from mpas_tools.io import open_dataset

from polaris.ocean.init_state import add_quiescent_normal_velocity
from polaris.ocean.model import OceanIOStep
from polaris.ocean.vertical import init_vertical_coord


class Init(OceanIOStep):
    """
    A step for building an idealized, quiescent (summer-stratified, zero
    velocity) initial condition on the combined Great Lakes mesh.

    No Great Lakes temperature/salinity climatology was available to build
    a realistic initial condition from, so the state is an idealized
    rest state instead (see docs/design_docs/great_lakes_mesh.md).

    Attributes
    ----------
    combine_step : polaris.tasks.ocean.great_lakes.combine.CombineGreatLakesMesh
        The step that builds the combined mesh this step initializes
    """  # noqa: E501

    def __init__(self, component, subdir, config, combine_step):
        """
        Create the step

        Parameters
        ----------
        component : polaris.tasks.ocean.Ocean
            The ocean component the step belongs to

        subdir : str
            The subdirectory for the step

        config : polaris.config.PolarisConfigParser
            A shared config parser

        combine_step : polaris.tasks.ocean.great_lakes.combine.CombineGreatLakesMesh
            The step that builds the combined mesh
        """  # noqa: E501
        super().__init__(
            component=component,
            name='init',
            subdir=subdir,
            ntasks=1,
            min_tasks=1,
        )
        self.set_shared_config(config, link='great_lakes_init.cfg')
        self.combine_step = combine_step

        self.add_input_file(
            filename='combined_mesh.nc',
            work_dir_target=f'{combine_step.path}/combined_mesh.nc',
        )
        self.add_input_file(
            filename='combined_topography.nc',
            work_dir_target=f'{combine_step.path}/combined_topography.nc',
        )
        self.add_input_file(
            filename='combined_graph.info',
            work_dir_target=f'{combine_step.path}/combined_graph.info',
        )
        self.add_input_file(
            filename='reconstruction_weights.nc',
            work_dir_target=(f'{combine_step.path}/reconstruction_weights.nc'),
        )

    def setup(self):
        """
        Declare output files based on the configured ocean model
        """
        super().setup()
        self.add_output_files_for_ocean_model_input(
            horiz_mesh_filename='mesh.nc',
            vert_coord_filename='vert_coord.nc',
            init_filename='init.nc',
            graph_filename='combined_graph.info',
        )

    def run(self):
        """
        Run this step
        """
        config = self.config
        section = config['great_lakes_init']
        min_bottom_depth = section.getfloat('min_bottom_depth')
        surface_temperature = section.getfloat('surface_temperature')
        bottom_temperature = section.getfloat('bottom_temperature')
        thermocline_depth = section.getfloat('thermocline_depth')
        thermocline_width = section.getfloat('thermocline_width')

        ds_mesh = open_dataset('combined_mesh.nc')
        ds_topo = open_dataset('combined_topography.nc')

        # one salinity value per basin (see great_lakes_init.cfg), looked
        # up by the basin_names/basinIndex that CombineGreatLakesMesh adds
        # to combined_topography.nc
        basin_names = ds_topo.attrs['basin_names'].split(',')
        salinity_by_basin = [
            section.getfloat(f'salinity_{name}') for name in basin_names
        ]

        # the merged mesh already carries fCell/fEdge/fVertex (added in
        # CombineGreatLakesMesh), unlike a bare culled mesh
        self.write_horiz_mesh_dataset(ds_mesh, 'mesh.nc', config)

        ds = ds_mesh.copy()
        ds['bottomDepth'] = _bottom_depth_from_elevation(
            ds_topo['base_elevation'], min_bottom_depth
        )
        ds['ssh'] = xr.zeros_like(ds.bottomDepth)
        init_vertical_coord(config, ds)

        # idealized summer-stratified rest state -- T/S must be set before
        # write_vert_coord_dataset, below: pseudothickness_from_ds()
        # silently returns None if temperature/salinity are absent, which
        # would leave RefPseudoThickness missing from vert_coord.nc
        ds['temperature'] = _idealized_temperature(
            ds.zMid.isel(Time=0),
            surface_temperature=surface_temperature,
            bottom_temperature=bottom_temperature,
            thermocline_depth=thermocline_depth,
            thermocline_width=thermocline_width,
        )
        ds['salinity'] = _salinity_by_basin(
            ds['temperature'], ds_topo['basinIndex'], salinity_by_basin
        )
        # mutates and returns ds -- not a DataArray to assign onto ds
        ds = add_quiescent_normal_velocity(ds, ds_mesh)

        self.write_vert_coord_dataset(ds, 'vert_coord.nc', config)
        self.write_initial_state_dataset(
            ds,
            'init.nc',
            config,
            lon=np.rad2deg(ds_mesh.lonCell.values),
            lat=np.rad2deg(ds_mesh.latCell.values),
        )


def _bottom_depth_from_elevation(base_elevation, min_bottom_depth):
    """
    Convert a negative-down bathymetric elevation to a positive-down
    ``bottomDepth``, clamping both fill-value/NaN cells and cells
    shallower than ``min_bottom_depth`` to that floor.

    Parameters
    ----------
    base_elevation : xarray.DataArray
        Bathymetric elevation, negative down, with NaN/fill-value cells
        where the remap had no valid source data

    min_bottom_depth : float
        The positive-down depth, in m, to clamp to

    Returns
    -------
    bottom_depth : xarray.DataArray
        The positive-down depth of the seafloor, with dimension ``nCells``
    """
    depth = -base_elevation
    depth = depth.where(np.isfinite(depth), min_bottom_depth)
    depth = np.maximum(depth, min_bottom_depth)
    depth.attrs = {
        'long_name': 'Depth of the bottom of the ocean',
        'units': 'm',
    }
    return depth.rename('bottomDepth')


def _salinity_by_basin(temperature, basin_index, salinity_by_basin):
    """
    A salinity field that is uniform with depth and within each basin, but
    varies between basins, since the Great Lakes' real background
    salinities differ enough between basins to be worth representing (see
    docs/design_docs/great_lakes_mesh.md).

    Parameters
    ----------
    temperature : xarray.DataArray
        The idealized temperature field, used only for its shape/dims so
        the returned salinity broadcasts against it directly

    basin_index : xarray.DataArray
        A per-cell index (dimension ``nCells``) into ``salinity_by_basin``

    salinity_by_basin : list of float
        The practical/absolute salinity, in g/kg, for each basin, ordered
        to match ``basin_index``

    Returns
    -------
    salinity : xarray.DataArray
        The salinity at each layer midpoint, with the same dimensions as
        ``temperature``
    """
    lookup = np.asarray(salinity_by_basin)
    salinity_per_cell = xr.DataArray(
        lookup[basin_index.values], dims=('nCells',)
    )
    salinity = xr.ones_like(temperature) * salinity_per_cell
    salinity.attrs = {
        'long_name': 'salinity',
        'units': 'g kg-1',
    }
    return salinity.rename('salinity')


def _idealized_temperature(
    z_mid,
    surface_temperature,
    bottom_temperature,
    thermocline_depth,
    thermocline_width,
):
    """
    An idealized summer-stratified temperature profile: uniform near-
    surface temperature, a tanh thermocline, and a cold, well-mixed
    bottom layer.

    Parameters
    ----------
    z_mid : xarray.DataArray
        The elevation (positive up) of layer midpoints, with dimensions
        ``(nCells, nVertLevels)``

    surface_temperature : float
        The temperature, in deg C, above the thermocline

    bottom_temperature : float
        The temperature, in deg C, below the thermocline

    thermocline_depth : float
        The positive-down depth, in m, of the thermocline's center

    thermocline_width : float
        The e-folding-like width, in m, of the tanh thermocline transition

    Returns
    -------
    temperature : xarray.DataArray
        The temperature at each layer midpoint, with a ``Time`` dimension
        added, ready to assign to an initial-condition dataset
    """
    depth = -z_mid
    frac_above_thermocline = 0.5 * (
        1.0 - np.tanh((depth - thermocline_depth) / thermocline_width)
    )
    temperature = (
        bottom_temperature
        + (surface_temperature - bottom_temperature) * frac_above_thermocline
    )
    temperature = temperature.transpose('nCells', 'nVertLevels')
    temperature = temperature.expand_dims(dim='Time', axis=0)
    temperature.attrs = {
        'long_name': 'potential temperature',
        'units': 'degrees Celsius',
    }
    return temperature.rename('temperature')
