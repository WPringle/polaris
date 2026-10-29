import itertools

import numpy as np
import xarray as xr
from mpas_tools.io import open_dataset, write_netcdf
from shapely import STRtree
from shapely.geometry import Polygon

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

        _check_no_basin_overlap(ds_meshes, lake_names)

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


def _check_no_basin_overlap(ds_meshes, lake_names):
    """
    Verify that no two basins' culled meshes physically overlap.

    Each basin is built independently and is assumed to be its own
    disjoint sub-domain (see
    :py:data:`polaris.tasks.ocean.great_lakes.lakes.GREAT_LAKES`); this
    checks that assumption instead of just asserting it, so a shared-
    geometry bug (e.g. two basins' bathymetry both reaching into the same
    connecting channel, as happened at the Soo Locks) fails loudly here
    rather than only being visible in a plot.

    Parameters
    ----------
    ds_meshes : list of xarray.Dataset
        Each basin's culled mesh, before ID offsetting

    lake_names : list of str
        The basin names, in the same order as ``ds_meshes``

    Raises
    ------
    ValueError
        If any two basins have cells that physically overlap
    """
    bounds = [_lon_lat_bounds(ds) for ds in ds_meshes]
    for i, j in itertools.combinations(range(len(ds_meshes)), 2):
        if not _bounds_intersect(bounds[i], bounds[j]):
            continue
        overlap = _overlapping_cells(
            ds_meshes[i], ds_meshes[j], bounds[i], bounds[j]
        )
        if overlap is not None:
            lon_min, lon_max, lat_min, lat_max = overlap['extent']
            raise ValueError(
                f'{lake_names[i]} and {lake_names[j]} have '
                f'{overlap["n_pairs"]} overlapping cell pairs '
                f'({overlap["n_i"]} {lake_names[i]} cells, '
                f'{overlap["n_j"]} {lake_names[j]} cells) spanning lon '
                f'[{lon_min:.4f}, {lon_max:.4f}], lat '
                f'[{lat_min:.4f}, {lat_max:.4f}]. Each basin is expected '
                f'to be its own disjoint sub-domain; add or adjust a '
                f'bathymetry clip/exclusion in [great_lakes_topo] to '
                f'separate them (see '
                f'docs/design_docs/great_lakes_mesh.md).'
            )


def _lon_lat_bounds(ds):
    """
    The (lon_min, lon_max, lat_min, lat_max) bounding box of a mesh's
    cell *polygons*, in degrees with longitude wrapped to ``[-180, 180)``.

    Uses vertex coordinates, not cell centers: a boundary cell's polygon
    can extend past its own center by roughly a cell width, which is
    exactly the scale at which two independently built basins can
    overlap, so a center-based bounding box could miss the very overlap
    this check exists to catch.
    """
    lon_deg = np.rad2deg(ds['lonVertex'].values)
    lon_deg = np.where(lon_deg > 180.0, lon_deg - 360.0, lon_deg)
    lat_deg = np.rad2deg(ds['latVertex'].values)
    return (
        float(lon_deg.min()),
        float(lon_deg.max()),
        float(lat_deg.min()),
        float(lat_deg.max()),
    )


def _bounds_intersect(bounds_a, bounds_b):
    """Whether two (lon_min, lon_max, lat_min, lat_max) boxes intersect"""
    a_lon_min, a_lon_max, a_lat_min, a_lat_max = bounds_a
    b_lon_min, b_lon_max, b_lat_min, b_lat_max = bounds_b
    return (
        a_lon_min <= b_lon_max
        and b_lon_min <= a_lon_max
        and a_lat_min <= b_lat_max
        and b_lat_min <= a_lat_max
    )


def _cell_polygons_in_box(ds, extent):
    """
    Build shapely polygons for the cells of ``ds`` whose own vertex
    bounding box intersects ``extent``.

    Selecting by each cell's own geometry, not its center, matters here:
    a boundary cell's polygon can extend past its center by roughly a
    cell width, which is exactly the scale of the overlap this check
    exists to catch, so a center-based filter could silently drop the
    one cell that actually overlaps.

    Parameters
    ----------
    ds : xarray.Dataset
        A basin's culled mesh

    extent : tuple of float
        The (lon_min, lon_max, lat_min, lat_max) box, in degrees

    Returns
    -------
    cell_indices : numpy.ndarray
        The indices of the selected cells

    polygons : list of shapely.geometry.Polygon
        One polygon per selected cell, in the same order as
        ``cell_indices``
    """
    lon_min, lon_max, lat_min, lat_max = extent
    lon_vertex = np.rad2deg(ds['lonVertex'].values)
    lon_vertex = np.where(lon_vertex > 180.0, lon_vertex - 360.0, lon_vertex)
    lat_vertex = np.rad2deg(ds['latVertex'].values)
    vertices_on_cell = ds['verticesOnCell'].values
    n_edges_on_cell = ds['nEdgesOnCell'].values
    max_edges = vertices_on_cell.shape[1]

    # verticesOnCell is one-based and zero-padded out to maxEdges, so the
    # padding is masked off (rather than indexed) before reducing
    valid = np.arange(max_edges)[None, :] < n_edges_on_cell[:, None]
    indices = np.where(valid, vertices_on_cell - 1, 0)
    cell_lon = lon_vertex[indices]
    cell_lat = lat_vertex[indices]
    cell_lon_min = np.where(valid, cell_lon, np.inf).min(axis=1)
    cell_lon_max = np.where(valid, cell_lon, -np.inf).max(axis=1)
    cell_lat_min = np.where(valid, cell_lat, np.inf).min(axis=1)
    cell_lat_max = np.where(valid, cell_lat, -np.inf).max(axis=1)

    mask = (
        (cell_lon_max >= lon_min)
        & (cell_lon_min <= lon_max)
        & (cell_lat_max >= lat_min)
        & (cell_lat_min <= lat_max)
    )
    cell_indices = np.flatnonzero(mask)

    polygons = []
    for cell_index in cell_indices:
        n = n_edges_on_cell[cell_index]
        vertex_indices = vertices_on_cell[cell_index, :n] - 1
        points = zip(
            lon_vertex[vertex_indices],
            lat_vertex[vertex_indices],
            strict=True,
        )
        polygons.append(Polygon(points))
    return cell_indices, polygons


def _overlapping_cells(ds_i, ds_j, bounds_i, bounds_j):
    """
    Find cells of two basins' meshes whose polygons physically overlap,
    restricted to the intersection of their bounding boxes.

    Parameters
    ----------
    ds_i, ds_j : xarray.Dataset
        Two basins' culled meshes

    bounds_i, bounds_j : tuple of float
        Each basin's (lon_min, lon_max, lat_min, lat_max) bounding box

    Returns
    -------
    overlap : dict or None
        ``None`` if no cells overlap, otherwise a dict with ``n_pairs``,
        ``n_i``, and ``n_j`` (numbers of overlapping cell pairs/cells) and
        ``extent`` (the (lon_min, lon_max, lat_min, lat_max) of the
        overlapping region)
    """
    extent = (
        max(bounds_i[0], bounds_j[0]),
        min(bounds_i[1], bounds_j[1]),
        max(bounds_i[2], bounds_j[2]),
        min(bounds_i[3], bounds_j[3]),
    )
    cell_indices_i, polygons_i = _cell_polygons_in_box(ds_i, extent)
    cell_indices_j, polygons_j = _cell_polygons_in_box(ds_j, extent)
    if len(polygons_i) == 0 or len(polygons_j) == 0:
        return None

    tree = STRtree(polygons_j)
    cells_i = set()
    cells_j = set()
    n_pairs = 0
    lon_min, lon_max, lat_min, lat_max = (1.0e9, -1.0e9, 1.0e9, -1.0e9)
    for local_i, polygon_i in enumerate(polygons_i):
        for local_j in tree.query(polygon_i, predicate='intersects'):
            overlap_area = polygon_i.intersection(polygons_j[local_j]).area
            if overlap_area > 0.0:
                n_pairs += 1
                cells_i.add(cell_indices_i[local_i])
                cells_j.add(cell_indices_j[local_j])
                bounds = polygon_i.intersection(polygons_j[local_j]).bounds
                lon_min = min(lon_min, bounds[0])
                lat_min = min(lat_min, bounds[1])
                lon_max = max(lon_max, bounds[2])
                lat_max = max(lat_max, bounds[3])

    if n_pairs == 0:
        return None
    return {
        'n_pairs': n_pairs,
        'n_i': len(cells_i),
        'n_j': len(cells_j),
        'extent': (lon_min, lon_max, lat_min, lat_max),
    }
