from polaris.tasks.ocean.great_lakes.lakes import GREAT_LAKES as GREAT_LAKES
from polaris.tasks.ocean.great_lakes.task import (
    GreatLakesMeshTask as GreatLakesMeshTask,
)


def add_great_lakes_tasks(component):
    """
    Add tasks for building a mesh for each Great Lake to the ocean
    component.

    Parameters
    ----------
    component : polaris.tasks.ocean.Ocean
        The ocean component to which the tasks will be added.
    """
    for lake_name in GREAT_LAKES:
        component.add_task(
            GreatLakesMeshTask(component=component, lake_name=lake_name)
        )
