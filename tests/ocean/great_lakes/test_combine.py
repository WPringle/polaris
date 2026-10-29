"""
Unit tests for :py:mod:`polaris.tasks.ocean.great_lakes.combine`.

All tests are self-contained: no file I/O and no full Polaris step
framework.  Small, hand-built MPAS mesh datasets stand in for two basins'
culled meshes, deliberately built with *different* ``maxEdges`` (as
Superior/Michigan-Huron's 9 and Erie/Ontario's 8 actually are), since that
mismatch is what ``mpas_tools.merge_grids.merge_grids`` cannot handle.
"""

import numpy as np
import pytest
import xarray as xr

from polaris.tasks.ocean.great_lakes.combine import (
    _apply_id_offsets,
    _check_no_basin_overlap,
    _merge_meshes,
    _pad_ragged_dims,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_basin_mesh(n_cells, n_edges, n_vertices, max_edges, n_real_edges=2):
    """
    A small, self-contained MPAS mesh dataset for one basin.

    Not a physically valid mesh (the point/connectivity geometry is not
    consistent), but every connectivity array only ever references other
    elements of the *same* basin (1-based, as real MPAS IDs are), with
    zero as the "no neighbor"/padding sentinel beyond ``n_real_edges`` --
    everything :py:mod:`polaris.tasks.ocean.great_lakes.combine` actually
    looks at.
    """
    max_edges2 = 2 * max_edges

    def _ring_ids(size, width, n_real, id_pool_size):
        # row i's "real" entries cycle through 1-based IDs from a pool of
        # size id_pool_size, padded with the zero sentinel out to width
        local = np.zeros((size, width), dtype=np.int32)
        idx = np.arange(size)
        for col in range(n_real):
            local[:, col] = ((idx + col) % id_pool_size) + 1
        return local

    ds = xr.Dataset()
    ds['indexToCellID'] = ('nCells', np.arange(1, n_cells + 1, dtype=np.int32))
    ds['indexToEdgeID'] = ('nEdges', np.arange(1, n_edges + 1, dtype=np.int32))
    ds['indexToVertexID'] = (
        'nVertices',
        np.arange(1, n_vertices + 1, dtype=np.int32),
    )
    ds['nEdgesOnCell'] = (
        'nCells',
        np.full(n_cells, n_real_edges, dtype=np.int32),
    )
    ds['cellsOnCell'] = (
        ('nCells', 'maxEdges'),
        _ring_ids(n_cells, max_edges, n_real_edges, n_cells),
    )
    ds['edgesOnCell'] = (
        ('nCells', 'maxEdges'),
        _ring_ids(n_cells, max_edges, n_real_edges, n_edges),
    )
    ds['verticesOnCell'] = (
        ('nCells', 'maxEdges'),
        _ring_ids(n_cells, max_edges, n_real_edges, n_vertices),
    )
    ds['cellsOnEdge'] = (
        ('nEdges', 'TWO'),
        _ring_ids(n_edges, 2, 2, n_cells),
    )
    ds['verticesOnEdge'] = (
        ('nEdges', 'TWO'),
        _ring_ids(n_edges, 2, 2, n_vertices),
    )
    edges_on_edge = _ring_ids(n_edges, max_edges2, n_real_edges, n_edges)
    ds['edgesOnEdge'] = (('nEdges', 'maxEdges2'), edges_on_edge)
    ds['weightsOnEdge'] = (
        ('nEdges', 'maxEdges2'),
        np.where(edges_on_edge > 0, 0.5, 0.0),
    )
    ds['edgesOnVertex'] = (
        ('nVertices', 'vertexDegree'),
        _ring_ids(n_vertices, 3, 3, n_edges),
    )
    ds['cellsOnVertex'] = (
        ('nVertices', 'vertexDegree'),
        _ring_ids(n_vertices, 3, 3, n_cells),
    )
    ds['latCell'] = ('nCells', np.linspace(0.0, 1.0, n_cells))
    ds['lonCell'] = ('nCells', np.linspace(0.0, 2.0, n_cells))
    ds.attrs.update(
        {
            'on_a_sphere': 'YES',
            'sphere_radius': 6371008.7714,
            'is_periodic': 'NO',
            'mesh_spec': '1.0',
            'Conventions': 'MPAS',
        }
    )
    return ds


# ---------------------------------------------------------------------------
# _pad_ragged_dims
# ---------------------------------------------------------------------------


def test_pad_ragged_dims_pads_with_zero():
    ds = _make_basin_mesh(
        n_cells=3, n_edges=3, n_vertices=3, max_edges=3, n_real_edges=2
    )
    padded = _pad_ragged_dims(ds, max_edges=5)

    assert padded.sizes['maxEdges'] == 5
    assert padded.sizes['maxEdges2'] == 10
    # the original (real) columns are untouched
    assert np.array_equal(
        padded['cellsOnCell'].values[:, :3], ds['cellsOnCell'].values
    )
    # the new columns are the zero sentinel, not garbage
    assert np.all(padded['cellsOnCell'].values[:, 3:] == 0)
    assert np.all(padded['edgesOnEdge'].values[:, 6:] == 0)
    assert np.all(padded['weightsOnEdge'].values[:, 6:] == 0.0)


def test_pad_ragged_dims_is_a_no_op_when_already_the_target_width():
    ds = _make_basin_mesh(n_cells=2, n_edges=2, n_vertices=2, max_edges=4)
    padded = _pad_ragged_dims(ds, max_edges=4)
    assert padded is ds


# ---------------------------------------------------------------------------
# _apply_id_offsets
# ---------------------------------------------------------------------------


def test_apply_id_offsets_shifts_nonzero_entries_only():
    ds = _make_basin_mesh(n_cells=3, n_edges=3, n_vertices=3, max_edges=3)
    offset = _apply_id_offsets(ds, cell_off=10, edge_off=20, vertex_off=30)

    assert np.array_equal(
        offset['indexToCellID'].values, np.array([11, 12, 13])
    )
    assert np.array_equal(
        offset['indexToEdgeID'].values, np.array([21, 22, 23])
    )
    assert np.array_equal(
        offset['indexToVertexID'].values, np.array([31, 32, 33])
    )

    # cell-ID connectivity shifted by cell_off, zero sentinel preserved
    coc = offset['cellsOnCell'].values
    assert np.all(coc[ds['cellsOnCell'].values == 0] == 0)
    assert np.array_equal(
        coc[ds['cellsOnCell'].values > 0],
        ds['cellsOnCell'].values[ds['cellsOnCell'].values > 0] + 10,
    )

    # edge-ID and vertex-ID connectivity shifted by their own offsets
    eoc = offset['edgesOnCell'].values
    assert np.array_equal(
        eoc[ds['edgesOnCell'].values > 0],
        ds['edgesOnCell'].values[ds['edgesOnCell'].values > 0] + 20,
    )
    voc = offset['verticesOnCell'].values
    assert np.array_equal(
        voc[ds['verticesOnCell'].values > 0],
        ds['verticesOnCell'].values[ds['verticesOnCell'].values > 0] + 30,
    )


# ---------------------------------------------------------------------------
# _merge_meshes
# ---------------------------------------------------------------------------


def test_merge_meshes_with_mismatched_max_edges():
    # Superior/Michigan-Huron use maxEdges=9, Erie/Ontario use 8 in the
    # real meshes; reproduce that mismatch here, since it is exactly what
    # mpas_tools.merge_grids.merge_grids cannot handle (a shape-mismatch
    # ValueError when concatenating two inputs with different maxEdges).
    ds_a = _make_basin_mesh(
        n_cells=4, n_edges=4, n_vertices=4, max_edges=4, n_real_edges=2
    )
    ds_b = _make_basin_mesh(
        n_cells=3, n_edges=3, n_vertices=3, max_edges=3, n_real_edges=2
    )

    ds_combined = _merge_meshes([ds_a, ds_b])

    assert ds_combined.sizes['nCells'] == 7
    assert ds_combined.sizes['nEdges'] == 7
    assert ds_combined.sizes['nVertices'] == 7
    assert ds_combined.sizes['maxEdges'] == 4
    assert ds_combined.sizes['maxEdges2'] == 8

    # indexToCellID/Edge/Vertex ID are a contiguous 1..N range
    assert np.array_equal(ds_combined['indexToCellID'].values, np.arange(1, 8))

    # no cross-basin connectivity: basin b's (second basin's) rows only
    # ever reference IDs in its own offset range, never basin a's
    coc = ds_combined['cellsOnCell'].values
    b_rows = coc[4:7]
    nonzero = b_rows[b_rows > 0]
    assert nonzero.min() >= 5
    assert nonzero.max() <= 7
    a_rows = coc[0:4]
    nonzero_a = a_rows[a_rows > 0]
    assert nonzero_a.min() >= 1
    assert nonzero_a.max() <= 4

    # basin b's padded columns (it only had maxEdges=3) are zero, not
    # garbage, after being widened to match basin a's maxEdges=4
    assert np.all(coc[4:7, 3:] == 0)

    # plain (non-connectivity) fields are concatenated, not altered
    assert np.allclose(
        ds_combined['latCell'].values,
        np.concatenate([ds_a['latCell'].values, ds_b['latCell'].values]),
    )

    # global attributes that agree between inputs are carried forward
    assert ds_combined.attrs['on_a_sphere'] == 'YES'
    assert ds_combined.attrs['Conventions'] == 'MPAS'


def test_merge_meshes_counts_equal_sum_of_inputs():
    ds_a = _make_basin_mesh(n_cells=5, n_edges=6, n_vertices=4, max_edges=6)
    ds_b = _make_basin_mesh(n_cells=2, n_edges=3, n_vertices=2, max_edges=6)
    ds_c = _make_basin_mesh(n_cells=3, n_edges=4, n_vertices=3, max_edges=6)

    ds_combined = _merge_meshes([ds_a, ds_b, ds_c])

    assert ds_combined.sizes['nCells'] == 10
    assert ds_combined.sizes['nEdges'] == 13
    assert ds_combined.sizes['nVertices'] == 9


# ---------------------------------------------------------------------------
# _check_no_basin_overlap
# ---------------------------------------------------------------------------


def _make_square_cell_mesh(lon_deg, lat_deg):
    """
    A one-cell basin mesh whose single cell is the quadrilateral with the
    given vertex (lon, lat) corners, in degrees, listed in order around
    the polygon.
    """
    lon_rad = np.deg2rad(np.asarray(lon_deg))
    lat_rad = np.deg2rad(np.asarray(lat_deg))
    ds = xr.Dataset()
    ds['lonVertex'] = ('nVertices', lon_rad)
    ds['latVertex'] = ('nVertices', lat_rad)
    ds['lonCell'] = ('nCells', np.array([lon_rad.mean()]))
    ds['latCell'] = ('nCells', np.array([lat_rad.mean()]))
    ds['nEdgesOnCell'] = ('nCells', np.array([4], dtype=np.int32))
    ds['verticesOnCell'] = (
        ('nCells', 'maxEdges'),
        np.array([[1, 2, 3, 4]], dtype=np.int32),
    )
    return ds


def test_check_no_basin_overlap_passes_for_disjoint_basins():
    ds_a = _make_square_cell_mesh([0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 1.0, 1.0])
    ds_b = _make_square_cell_mesh(
        [10.0, 11.0, 11.0, 10.0], [10.0, 10.0, 11.0, 11.0]
    )
    # should not raise
    _check_no_basin_overlap([ds_a, ds_b], ['a', 'b'])


def test_check_no_basin_overlap_passes_for_touching_basins():
    # sharing an edge (touching, zero-area intersection) is not an overlap
    ds_a = _make_square_cell_mesh([0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 1.0, 1.0])
    ds_b = _make_square_cell_mesh([1.0, 2.0, 2.0, 1.0], [0.0, 0.0, 1.0, 1.0])
    _check_no_basin_overlap([ds_a, ds_b], ['a', 'b'])


def test_check_no_basin_overlap_raises_for_overlapping_basins():
    ds_a = _make_square_cell_mesh([0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 1.0, 1.0])
    ds_b = _make_square_cell_mesh([0.5, 1.5, 1.5, 0.5], [0.5, 0.5, 1.5, 1.5])
    with pytest.raises(ValueError, match='overlapping cell pairs'):
        _check_no_basin_overlap([ds_a, ds_b], ['a', 'b'])
