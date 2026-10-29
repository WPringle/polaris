"""
Unit tests for :py:mod:`polaris.tasks.ocean.great_lakes.viz`.

All tests are self-contained: no file I/O and no full Polaris step
framework.
"""

import numpy as np
import xarray as xr

from polaris.tasks.ocean.great_lakes.viz import _cell_width_km, _cells_in_box


def _make_mesh(lon_deg, lat_deg, dc_edge_km, n_edges_on_cell, edges_on_cell):
    """A minimal MPAS mesh dataset with just the fields these helpers use"""
    ds = xr.Dataset()
    ds['lonCell'] = ('nCells', np.deg2rad(np.asarray(lon_deg)))
    ds['latCell'] = ('nCells', np.deg2rad(np.asarray(lat_deg)))
    ds['dcEdge'] = ('nEdges', 1.0e3 * np.asarray(dc_edge_km))
    ds['nEdgesOnCell'] = (
        'nCells',
        np.asarray(n_edges_on_cell, dtype=np.int32),
    )
    ds['edgesOnCell'] = (
        ('nCells', 'maxEdges'),
        np.asarray(edges_on_cell, dtype=np.int32),
    )
    return ds


def test_cell_width_km_masks_padding():
    # 3 cells with ragged nEdgesOnCell (5, 6, 3) and maxEdges = 6; the
    # zero-padded slots of edgesOnCell must not contribute to the mean
    dc_edge_km = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    edges_on_cell = [
        [1, 2, 3, 4, 5, 0],
        [1, 2, 3, 4, 5, 6],
        [1, 2, 3, 0, 0, 0],
    ]
    n_edges_on_cell = [5, 6, 3]
    ds_mesh = _make_mesh(
        lon_deg=[0.0, 0.0, 0.0],
        lat_deg=[0.0, 0.0, 0.0],
        dc_edge_km=dc_edge_km,
        n_edges_on_cell=n_edges_on_cell,
        edges_on_cell=edges_on_cell,
    )

    width_km = _cell_width_km(ds_mesh)

    np.testing.assert_allclose(width_km[0], np.mean(dc_edge_km[:5]))
    np.testing.assert_allclose(width_km[1], np.mean(dc_edge_km[:6]))
    np.testing.assert_allclose(width_km[2], np.mean(dc_edge_km[:3]))


def test_cell_width_km_uniform_dc_edge():
    # a cell whose edges all share one dcEdge returns exactly that value,
    # also checking the m -> km conversion
    edges_on_cell = [[1, 2, 3, 4]]
    ds_mesh = _make_mesh(
        lon_deg=[0.0],
        lat_deg=[0.0],
        dc_edge_km=[2.5, 2.5, 2.5, 2.5],
        n_edges_on_cell=[4],
        edges_on_cell=edges_on_cell,
    )

    width_km = _cell_width_km(ds_mesh)

    np.testing.assert_allclose(width_km, [2.5])


def test_cells_in_box():
    # one cell inside the box, one just outside each of the four faces,
    # and one whose longitude is stored past the +/-180 wrap
    lon_deg = [-81.0, -92.5, -75.7, -81.0, -81.0, 170.0]
    lat_deg = [45.0, 45.0, 45.0, 41.1, 49.8, 45.0]
    ds_mesh = _make_mesh(
        lon_deg=[lon % 360.0 for lon in lon_deg],
        lat_deg=lat_deg,
        dc_edge_km=[1.0],
        n_edges_on_cell=[1] * len(lon_deg),
        edges_on_cell=[[1]] * len(lon_deg),
    )
    extent = (-92.4, -75.8, 41.2, 49.7)

    cell_indices = _cells_in_box(ds_mesh, extent)

    np.testing.assert_array_equal(cell_indices, [0])
