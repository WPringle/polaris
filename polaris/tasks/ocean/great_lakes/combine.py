import numpy as np
import xarray as xr
from mpas_tools.io import open_dataset, write_netcdf

from polaris import Step
from polaris.coriolis import add_coriolis_to_dataset
from polaris.mesh.reconstruct import (
    compute_reconstruction_weights,
    get_reconstruction_validate_vars,
)
from polaris.mesh.scrip import write_mpas_scrip_file
from polaris.mesh.validate import MPAS_MESH_VALIDATE_VARS
from polaris.model_step import make_graph_file

#: Global attributes carried forward from the first basin's mesh, since
#: all four basins are built from the same base-mesh generator and these
#: are expected to agree.
_GLOBAL_ATTRS = (
    'on_a_sphere',
    'sphere_radius',
    'is_periodic',
    'mesh_spec',
    'Conventions',
)

#: Connectivity variables whose entries are cell indices and must be
#: shifted by each basin's running cell-ID offset.  A zero entry is the
#: MPAS "no neighbor" sentinel and is left untouched.
_CELL_ID_VARS = ('cellsOnCell', 'cellsOnEdge', 'cellsOnVertex')

#: Connectivity variables whose entries are edge indices.
_EDGE_ID_VARS = ('edgesOnCell', 'edgesOnEdge', 'edgesOnVertex')

#: Connectivity variables whose entries are vertex indices.
_VERTEX_ID_VARS = ('verticesOnCell', 'verticesOnEdge')

#: Connectivity variables padded with zeros (the same "no neighbor"/
#: "no edge" sentinel used elsewhere in the ragged dimension) out to the
#: widest basin's ``maxEdges`` before concatenation.
_MAX_EDGES_VARS = ('cellsOnCell', 'edgesOnCell', 'verticesOnCell')

#: Connectivity/weight variables padded out to the widest basin's
#: ``maxEdges2`` (always ``2 * maxEdges``) before concatenation.
_MAX_EDGES2_VARS = ('edgesOnEdge', 'weightsOnEdge')


class CombineGreatLakesMesh(Step):
    """
    A step for merging each Great Lakes basin's culled mesh and remapped
    bathymetry into a single mesh spanning all basins.

    Each basin is already its own closed sub-domain (see
    :py:data:`polaris.tasks.ocean.great_lakes.lakes.GREAT_LAKES`), so the
    merge is a disjoint union: the basins' cell/edge/vertex arrays are
    concatenated with per-basin ID offsets and no connectivity is created
    between them.  This step also adds the Coriolis parameter, vector-
    reconstruction weights, and a SCRIP file for the merged mesh, since a
    spherical mesh needs all three before it can be used as an ocean
    model input (see
    :py:meth:`polaris.tasks.ocean.Ocean.write_horiz_mesh_dataset`).

    Attributes
    ----------
    lake_steps : dict of {str: dict}
        For each basin name, a dict with the basin's ``cull_mesh``
        (:py:class:`~polaris.tasks.ocean.great_lakes.cull.GreatLakesCullMesh`)
        and ``topo_remap``
        (:py:class:`~polaris.tasks.ocean.great_lakes.topo.remap.GreatLakesTopoRemap`)
        steps
    """

    def __init__(self, component, subdir, config, lake_steps):
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

        lake_steps : dict of {str: dict}
            For each basin name, a dict with the basin's ``cull_mesh`` and
            ``topo_remap`` steps
        """
        super().__init__(component=component, name='combine', subdir=subdir)
        self.set_shared_config(config, link='great_lakes_combined.cfg')
        self.lake_steps = lake_steps

        for lake_name, steps in lake_steps.items():
            self.add_input_file(
                filename=f'{lake_name}_mesh.nc',
                work_dir_target=f'{steps["cull_mesh"].path}/culled_mesh.nc',
            )
            self.add_input_file(
                filename=f'{lake_name}_topo.nc',
                work_dir_target=(
                    f'{steps["topo_remap"].path}/topography_remapped.nc'
                ),
            )

        self.add_output_file(
            filename='combined_mesh.nc', validate_vars=MPAS_MESH_VALIDATE_VARS
        )
        self.add_output_file(filename='combined_topography.nc')
        self.add_output_file(filename='combined_graph.info')
        self.add_output_file(
            filename='reconstruction_weights.nc',
            validate_vars=get_reconstruction_validate_vars(location='cell'),
        )
        self.add_output_file(filename='combined_mesh.scrip.nc')

    def run(self):
        """
        Run this step
        """
        config = self.config
        lake_names = list(self.lake_steps.keys())

        ds_meshes = [
            open_dataset(f'{lake_name}_mesh.nc') for lake_name in lake_names
        ]
        ds_topos = [
            open_dataset(f'{lake_name}_topo.nc') for lake_name in lake_names
        ]

        ds_combined = _merge_meshes(ds_meshes)
        ds_combined = add_coriolis_to_dataset(config, ds_combined)
        write_netcdf(ds_combined, 'combined_mesh.nc')

        ds_topo_combined = xr.concat(ds_topos, dim='nCells')
        ds_topo_combined['basinIndex'] = (
            'nCells',
            np.concatenate(
                [
                    np.full(ds.sizes['nCells'], index, dtype=np.int32)
                    for index, ds in enumerate(ds_topos)
                ]
            ),
        )
        ds_topo_combined['basinIndex'].attrs['long_name'] = (
            'index into the comma-separated basin_names global attribute '
            'identifying which basin each cell belongs to'
        )
        ds_topo_combined.attrs['basin_names'] = ','.join(lake_names)
        write_netcdf(ds_topo_combined, 'combined_topography.nc')

        make_graph_file(
            mesh_filename='combined_mesh.nc',
            graph_filename='combined_graph.info',
        )

        self.logger.info(
            'Compute vector-reconstruction weights at cell centers for '
            'the combined mesh'
        )
        ds_weights = compute_reconstruction_weights(
            ds_combined, location='cell'
        )
        write_netcdf(ds_weights, 'reconstruction_weights.nc')

        self.logger.info('Create SCRIP file for the combined mesh')
        write_mpas_scrip_file(
            mesh_filename='combined_mesh.nc',
            scrip_filename='combined_mesh.scrip.nc',
            mesh_name='great_lakes_combined',
            logger=self.logger,
        )


def _merge_meshes(ds_meshes):
    """
    Merge a list of culled MPAS meshes, each already its own disjoint
    closed sub-domain, into one mesh with renumbered global IDs and no
    connectivity between the inputs.

    Parameters
    ----------
    ds_meshes : list of xarray.Dataset
        The basins' culled meshes, in the order their cells should appear
        in the combined mesh

    Returns
    -------
    ds_combined : xarray.Dataset
        The combined mesh
    """
    max_edges = max(ds.sizes['maxEdges'] for ds in ds_meshes)
    ds_meshes = [_pad_ragged_dims(ds, max_edges) for ds in ds_meshes]

    n_cells = [ds.sizes['nCells'] for ds in ds_meshes]
    n_edges = [ds.sizes['nEdges'] for ds in ds_meshes]
    n_vertices = [ds.sizes['nVertices'] for ds in ds_meshes]
    cell_offsets = np.cumsum([0, *n_cells[:-1]])
    edge_offsets = np.cumsum([0, *n_edges[:-1]])
    vertex_offsets = np.cumsum([0, *n_vertices[:-1]])

    ds_offset = []
    for ds, cell_off, edge_off, vertex_off in zip(
        ds_meshes, cell_offsets, edge_offsets, vertex_offsets, strict=True
    ):
        ds_offset.append(_apply_id_offsets(ds, cell_off, edge_off, vertex_off))

    ds0 = ds_offset[0]
    cell_vars = [v for v in ds0.data_vars if 'nCells' in ds0[v].dims]
    edge_vars = [v for v in ds0.data_vars if 'nEdges' in ds0[v].dims]
    vertex_vars = [v for v in ds0.data_vars if 'nVertices' in ds0[v].dims]

    ds_cells = xr.concat([ds[cell_vars] for ds in ds_offset], dim='nCells')
    ds_edges = xr.concat([ds[edge_vars] for ds in ds_offset], dim='nEdges')
    ds_vertices = xr.concat(
        [ds[vertex_vars] for ds in ds_offset], dim='nVertices'
    )
    ds_combined = xr.merge(
        [ds_cells, ds_edges, ds_vertices], combine_attrs='override'
    )

    for attr in _GLOBAL_ATTRS:
        values = {ds.attrs.get(attr) for ds in ds_meshes}
        if len(values) == 1:
            ds_combined.attrs[attr] = values.pop()

    return ds_combined


def _pad_ragged_dims(ds, max_edges):
    """
    Pad a mesh's ragged ``maxEdges``/``maxEdges2`` connectivity variables
    with zeros (the existing MPAS "no neighbor"/"no edge" sentinel) out to
    a shared width, so that basins built with different ``maxEdges`` (as
    Superior/Michigan-Huron's 9 versus Erie/Ontario's 8 currently are) can
    be concatenated.  ``mpas_tools.merge_grids.merge_grids`` does not
    handle this: it pads the output dimension to the max of its two
    inputs but then assigns each input's un-padded array directly into
    the padded output, which raises a shape-mismatch error whenever the
    two inputs' ``maxEdges`` actually differ.
    """
    if ds.sizes['maxEdges'] == max_edges:
        return ds

    updates = {}
    for var in _MAX_EDGES_VARS:
        updates[var] = _pad_last_axis(ds[var], max_edges)
    for var in _MAX_EDGES2_VARS:
        updates[var] = _pad_last_axis(ds[var], 2 * max_edges)
    return ds.assign(updates)


def _pad_last_axis(da, target_size):
    """Pad a 2D data array's last axis with zeros out to ``target_size``"""
    pad = target_size - da.sizes[da.dims[-1]]
    values = np.pad(da.values, ((0, 0), (0, pad)), constant_values=0)
    return (da.dims, values.astype(da.dtype))


def _apply_id_offsets(ds, cell_off, edge_off, vertex_off):
    """
    Shift one basin's cell/edge/vertex IDs (and every connectivity array
    that references them) by its running offset in the combined mesh,
    following the same per-variable offset table as
    ``mpas_tools.merge_grids.merge_grids`` (zero entries are the MPAS
    "no neighbor" sentinel and are left untouched).
    """
    updates = {
        'indexToCellID': ds['indexToCellID'] + cell_off,
        'indexToEdgeID': ds['indexToEdgeID'] + edge_off,
        'indexToVertexID': ds['indexToVertexID'] + vertex_off,
    }
    for var, offset in [
        *((v, cell_off) for v in _CELL_ID_VARS),
        *((v, edge_off) for v in _EDGE_ID_VARS),
        *((v, vertex_off) for v in _VERTEX_ID_VARS),
    ]:
        values = ds[var].values.copy()
        mask = values > 0
        values[mask] += offset
        updates[var] = (ds[var].dims, values)
    return ds.assign(updates)
