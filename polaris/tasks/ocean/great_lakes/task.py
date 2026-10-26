import os

from polaris.config import PolarisConfigParser
from polaris.task import Task
from polaris.tasks.ocean.great_lakes.cull import GreatLakesCullMesh
from polaris.tasks.ocean.great_lakes.mesh import GreatLakesSphericalMesh
from polaris.tasks.ocean.great_lakes.topo import (
    CoarsenBathymetryStep,
    GreatLakesTopoMap,
    GreatLakesTopoRemap,
    PrepareBathymetryStep,
    WetMaskOnMeshStep,
)


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
