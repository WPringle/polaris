import os

from polaris.config import PolarisConfigParser
from polaris.task import Task
from polaris.tasks.ocean.great_lakes.combine import CombineGreatLakesMesh
from polaris.tasks.ocean.great_lakes.cull import GreatLakesCullMesh
from polaris.tasks.ocean.great_lakes.forward import Forward
from polaris.tasks.ocean.great_lakes.init import Init
from polaris.tasks.ocean.great_lakes.mesh import GreatLakesSphericalMesh
from polaris.tasks.ocean.great_lakes.topo import (
    CoarsenBathymetryStep,
    GreatLakesTopoMap,
    GreatLakesTopoRemap,
    PrepareBathymetryStep,
    WetMaskOnMeshStep,
)

# Lake Superior's downloaded bathymetry extends past the Soo Locks at
# Sault Ste. Marie into the St. Marys River approach/Whitefish Bay, water
# whose level is hydraulically controlled by the locks rather than part of
# the open lake -- see docs/design_docs/great_lakes_mesh.md.
_LON_MAX_DEG_OVERRIDES = {'superior': -84.34}


class GreatLakesMeshTask(Task):
    """
    A task for building a culled, bathymetry-informed mesh for a Great
    Lakes basin, which may be a single lake or (for ``michigan_huron``) a
    combined domain built from more than one downloaded bathymetry file
    (see :py:class:`polaris.tasks.ocean.great_lakes.topo.bathymetry.PrepareBathymetryStep`).

    The remaining connecting channels are not covered by any of the
    downloaded bathymetry, so each basin is otherwise its own closed
    sub-domain for now (see
    :py:data:`polaris.tasks.ocean.great_lakes.lakes.GREAT_LAKES`).

    Attributes
    ----------
    cull_mesh : polaris.tasks.ocean.great_lakes.cull.GreatLakesCullMesh
        The step that culls the basin's mesh down to its wet cells.
        Exposed so that
        :py:class:`polaris.tasks.ocean.great_lakes.combine.CombineGreatLakesMesh`
        can depend on it across basins without re-creating it.

    topo_remap_culled : polaris.tasks.ocean.great_lakes.topo.remap.GreatLakesTopoRemap
        The step that remaps the basin's bathymetry onto its culled mesh.
        Exposed for the same reason as ``cull_mesh``.
    """  # noqa: E501

    def __init__(self, component, lake_name):
        """
        Create the task

        Parameters
        ----------
        component : polaris.tasks.ocean.Ocean
            The ocean component the task belongs to

        lake_name : str
            The name of the lake
        """
        subdir = os.path.join('spherical', 'great_lakes', lake_name, 'mesh')
        super().__init__(
            component=component,
            name=f'great_lakes_{lake_name}_mesh',
            subdir=subdir,
        )

        config = PolarisConfigParser(
            filepath=os.path.join(
                component.name, subdir, 'great_lakes_mesh.cfg'
            )
        )
        config.add_from_package('polaris.mesh.spherical', 'spherical.cfg')
        config.add_from_package('polaris.remap', 'mapping.cfg')
        config.add_from_package(
            'polaris.tasks.ocean.great_lakes', 'great_lakes_mesh.cfg'
        )
        config.add_from_package(
            'polaris.tasks.ocean.great_lakes', 'great_lakes_topo.cfg'
        )
        config.set('spherical_mesh', 'mpas_mesh_filename', 'base_mesh.nc')
        if lake_name in _LON_MAX_DEG_OVERRIDES:
            config.set(
                'great_lakes_topo',
                'lon_max_deg',
                str(_LON_MAX_DEG_OVERRIDES[lake_name]),
            )
        self.set_shared_config(config, link='great_lakes_mesh.cfg')

        bathymetry_step = PrepareBathymetryStep(
            component=component,
            subdir=os.path.join(subdir, 'bathymetry'),
            config=config,
            lake_name=lake_name,
        )
        self.add_step(bathymetry_step)

        base_mesh = GreatLakesSphericalMesh(
            component=component,
            subdir=os.path.join(subdir, 'base_mesh'),
            bathymetry_step=bathymetry_step,
            lake_name=lake_name,
        )
        base_mesh.set_shared_config(config, link='great_lakes_mesh.cfg')
        self.add_step(base_mesh)

        wet_mask_step = WetMaskOnMeshStep(
            component=component,
            subdir=os.path.join(subdir, 'wet_mask'),
            config=config,
            bathymetry_step=bathymetry_step,
            base_mesh=base_mesh,
        )
        self.add_step(wet_mask_step)

        cull_mesh = GreatLakesCullMesh(
            component=component,
            subdir=os.path.join(subdir, 'cull_mesh'),
            config=config,
            base_mesh=base_mesh,
            wet_mask_step=wet_mask_step,
        )
        self.add_step(cull_mesh)
        self.cull_mesh = cull_mesh

        coarsen_step = CoarsenBathymetryStep(
            component=component,
            subdir=os.path.join(subdir, 'topo', 'coarsen'),
            config=config,
            bathymetry_step=bathymetry_step,
            lake_name=lake_name,
        )
        self.add_step(coarsen_step)

        topo_map_culled = GreatLakesTopoMap(
            component=component,
            name='topo_map_culled',
            subdir=os.path.join(subdir, 'topo', 'map_culled'),
            config=config,
            mesh_name=lake_name,
            bathymetry_step=coarsen_step,
            mesh_step=cull_mesh,
            mesh_filename='culled_mesh.nc',
            bathymetry_filename='bathymetry_coarse.nc',
        )
        self.add_step(topo_map_culled)

        topo_remap_culled = GreatLakesTopoRemap(
            component=component,
            name='topo_remap_culled',
            subdir=os.path.join(subdir, 'topo', 'remap_culled'),
            config=config,
            bathymetry_step=coarsen_step,
            topo_map=topo_map_culled,
            bathymetry_filename='bathymetry_coarse.nc',
        )
        self.add_step(topo_remap_culled)
        self.topo_remap_culled = topo_remap_culled


class GreatLakesCombinedMeshTask(Task):
    """
    A task that merges all four Great Lakes basins' culled meshes (and
    their remapped bathymetry) into a single mesh spanning all basins, and
    adds the Coriolis parameter, vector-reconstruction weights, and a
    SCRIP file the merged mesh needs before it can be used as an ocean
    model input.  See
    :py:class:`polaris.tasks.ocean.great_lakes.combine.CombineGreatLakesMesh`
    for why this is a disjoint union rather than a connected mesh.

    Attributes
    ----------
    combine_step : polaris.tasks.ocean.great_lakes.combine.CombineGreatLakesMesh
        The step that builds the combined mesh. Exposed so that the
        combined ``init`` task can depend on it.
    """  # noqa: E501

    def __init__(self, component, lake_tasks):
        """
        Create the task

        Parameters
        ----------
        component : polaris.tasks.ocean.Ocean
            The ocean component the task belongs to

        lake_tasks : dict of {str: GreatLakesMeshTask}
            Each basin's mesh task, keyed by basin name, already added to
            ``component``
        """
        subdir = os.path.join('spherical', 'great_lakes', 'combined', 'mesh')
        super().__init__(
            component=component,
            name='great_lakes_combined_mesh',
            subdir=subdir,
        )

        config = PolarisConfigParser(
            filepath=os.path.join(
                component.name, subdir, 'great_lakes_combined.cfg'
            )
        )
        config.add_from_package(
            'polaris.tasks.ocean.great_lakes', 'great_lakes_combined.cfg'
        )
        self.set_shared_config(config, link='great_lakes_combined.cfg')

        lake_steps = {
            lake_name: {
                'cull_mesh': task.cull_mesh,
                'topo_remap': task.topo_remap_culled,
            }
            for lake_name, task in lake_tasks.items()
        }

        combine_step = CombineGreatLakesMesh(
            component=component,
            subdir=os.path.join(subdir, 'combine'),
            config=config,
            lake_steps=lake_steps,
        )
        self.add_step(combine_step)
        self.combine_step = combine_step


class GreatLakesCombinedInitTask(Task):
    """
    A task that builds an idealized, quiescent initial condition on the
    combined Great Lakes mesh. See
    :py:class:`polaris.tasks.ocean.great_lakes.init.Init` for why the
    state is idealized rather than a realistic climatology.

    Attributes
    ----------
    init_step : polaris.tasks.ocean.great_lakes.init.Init
        The step that builds the initial condition. Exposed so that the
        combined ``forward`` task can depend on it.
    """

    def __init__(self, component, mesh_task):
        """
        Create the task

        Parameters
        ----------
        component : polaris.tasks.ocean.Ocean
            The ocean component the task belongs to

        mesh_task : GreatLakesCombinedMeshTask
            The task that builds the combined mesh, already added to
            ``component``
        """
        subdir = os.path.join('spherical', 'great_lakes', 'combined', 'init')
        super().__init__(
            component=component,
            name='great_lakes_combined_init',
            subdir=subdir,
        )

        config = PolarisConfigParser(
            filepath=os.path.join(
                component.name, subdir, 'great_lakes_init.cfg'
            )
        )
        config.add_from_package(
            'polaris.tasks.ocean.great_lakes', 'great_lakes_init.cfg'
        )
        self.set_shared_config(config, link='great_lakes_init.cfg')

        init_step = Init(
            component=component,
            subdir=os.path.join(subdir, 'init'),
            config=config,
            combine_step=mesh_task.combine_step,
        )
        self.add_step(init_step)
        self.init_step = init_step


class GreatLakesCombinedForwardTask(Task):
    """
    A task that runs a short forward smoke test on the combined Great
    Lakes mesh and idealized initial condition. See
    :py:class:`polaris.tasks.ocean.great_lakes.forward.Forward` for why
    this validates the pipeline rather than producing a tuned simulation.
    """

    def __init__(self, component, init_task):
        """
        Create the task

        Parameters
        ----------
        component : polaris.tasks.ocean.Ocean
            The ocean component the task belongs to

        init_task : GreatLakesCombinedInitTask
            The task that builds the initial condition, already added to
            ``component``
        """
        subdir = os.path.join(
            'spherical', 'great_lakes', 'combined', 'forward'
        )
        super().__init__(
            component=component,
            name='great_lakes_combined_forward',
            subdir=subdir,
        )

        config = PolarisConfigParser(
            filepath=os.path.join(
                component.name, subdir, 'great_lakes_forward.cfg'
            )
        )
        # for [ocean] eos_type: the forward run's EOS must match the one
        # the initial condition's tracers were written with
        config.add_from_package(
            'polaris.tasks.ocean.great_lakes', 'great_lakes_init.cfg'
        )
        config.add_from_package(
            'polaris.tasks.ocean.great_lakes', 'great_lakes_forward.cfg'
        )
        self.set_shared_config(config, link='great_lakes_forward.cfg')

        forward_step = Forward(
            component=component,
            subdir=os.path.join(subdir, 'forward'),
            config=config,
            init=init_task.init_step,
        )
        self.add_step(forward_step)
