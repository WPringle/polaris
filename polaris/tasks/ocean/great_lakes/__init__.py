from polaris.tasks.ocean.great_lakes.lakes import GREAT_LAKES as GREAT_LAKES
from polaris.tasks.ocean.great_lakes.task import (
    GreatLakesCombinedForwardTask as GreatLakesCombinedForwardTask,
)
from polaris.tasks.ocean.great_lakes.task import (
    GreatLakesCombinedInitTask as GreatLakesCombinedInitTask,
)
from polaris.tasks.ocean.great_lakes.task import (
    GreatLakesCombinedMeshTask as GreatLakesCombinedMeshTask,
)
from polaris.tasks.ocean.great_lakes.task import (
    GreatLakesMeshTask as GreatLakesMeshTask,
)


def add_great_lakes_tasks(component):
    """
    Add tasks for building a mesh for each Great Lake, for merging all
    four into a single combined mesh, and for building and running an
    idealized initial condition on that combined mesh, to the ocean
    component.

    Parameters
    ----------
    component : polaris.tasks.ocean.Ocean
        The ocean component to which the tasks will be added.
    """
    lake_tasks = {}
    for lake_name in GREAT_LAKES:
        task = GreatLakesMeshTask(component=component, lake_name=lake_name)
        component.add_task(task)
        lake_tasks[lake_name] = task

    mesh_task = GreatLakesCombinedMeshTask(
        component=component, lake_tasks=lake_tasks
    )
    component.add_task(mesh_task)

    init_task = GreatLakesCombinedInitTask(
        component=component, mesh_task=mesh_task
    )
    component.add_task(init_task)

    component.add_task(
        GreatLakesCombinedForwardTask(component=component, init_task=init_task)
    )
