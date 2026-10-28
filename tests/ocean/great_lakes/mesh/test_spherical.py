import numpy as np

from polaris.tasks.ocean.great_lakes.mesh.spherical import (
    _fill_enclosed_region,
    _length_scale_cell_width,
    _topographic_length_scale,
)

_DX_M = 100.0
_DY_M = 100.0
_THRESHOLD_HIGH = 0.1
_BOX = (16, 25, 16, 25)
_CENTER = (20, 20)


def _box_ring_mask(shape, box, gap=None, gap_width=1):
    """
    A boolean mask that is ``True`` along the perimeter of ``box`` (a
    (row_start, row_stop, col_start, col_stop) tuple) and ``False``
    everywhere else, with ``gap_width`` cells of the top edge cleared
    starting at column ``gap`` if given, to test that closing can bridge a
    small gap in an otherwise-enclosing wall.
    """
    mask = np.zeros(shape, dtype=bool)
    row_start, row_stop, col_start, col_stop = box
    mask[row_start, col_start:col_stop] = True
    mask[row_stop - 1, col_start:col_stop] = True
    mask[row_start:row_stop, col_start] = True
    mask[row_start:row_stop, col_stop - 1] = True
    if gap is not None:
        mask[row_start, gap : gap + gap_width] = False
    return mask


def test_fill_enclosed_region_fills_a_fully_enclosed_pocket():
    mask = _box_ring_mask(shape=(41, 41), box=_BOX)

    # too small a radius leaves the interior pocket unfilled
    closed_small = _fill_enclosed_region(
        mask=mask, dx_m=_DX_M, dy_m=_DY_M, radius_km=0.2
    )
    assert not closed_small[_CENTER]

    # a large-enough radius fills the interior, fully ringed by True
    # cells, without spilling out far beyond the ring
    closed_large = _fill_enclosed_region(
        mask=mask, dx_m=_DX_M, dy_m=_DY_M, radius_km=0.4
    )
    assert closed_large[_CENTER]
    assert not closed_large[0, 0]
    assert not closed_large[40, 40]


def test_fill_enclosed_region_bridges_a_gap_in_the_wall():
    mask = _box_ring_mask(shape=(41, 41), box=_BOX, gap=19, gap_width=3)
    closed = _fill_enclosed_region(
        mask=mask, dx_m=_DX_M, dy_m=_DY_M, radius_km=0.4
    )
    assert closed[_CENTER]


def test_fill_enclosed_region_does_not_fill_an_open_line():
    # a straight, open line (e.g. an ordinary shoreline drop-off, or a
    # long unenclosed stretch of short topographic length scale) encloses
    # nothing, so closing should not balloon out to fill either side of it
    mask = np.zeros((41, 41), dtype=bool)
    mask[20, :] = True
    closed = _fill_enclosed_region(
        mask=mask, dx_m=_DX_M, dy_m=_DY_M, radius_km=0.4
    )
    assert not closed[0, 20]
    assert not closed[40, 20]


def test_topographic_length_scale_flat_region_is_infinite():
    # a uniform-depth interior has no local gradient anywhere
    elevation = np.full((3, 5), -50.0)
    ocean_mask = np.full((3, 5), True)
    length_scale = _topographic_length_scale(
        elevation, ocean_mask, _DX_M, _DY_M
    )
    assert np.all(np.isinf(length_scale))


def test_topographic_length_scale_is_zero_next_to_land():
    # he = h1*h2/(h1+h2) collapses to zero whenever either neighbor is
    # land, regardless of how gentle the real underwater slope is -- this
    # is the near-coast artifact that production code excludes separately
    # via the coastal-transition buffer
    elevation = np.array(
        [
            [1.0, 1.0, 1.0, 1.0],
            [1.0, -10.0, -10.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
        ]
    )
    ocean_mask = elevation < 0.0
    length_scale = _topographic_length_scale(
        elevation, ocean_mask, _DX_M, _DY_M
    )
    assert length_scale[1, 1] == 0.0
    assert length_scale[1, 2] == 0.0


def test_topographic_length_scale_matches_pairwise_formula():
    # a symmetric depth outlier (10 m) surrounded by flat 50 m depth on
    # both sides, far from any land, isolates a single pairwise
    # computation at the outlier and its two neighbors
    depths = [-50.0, -50.0, -50.0, -10.0, -50.0, -50.0, -50.0]
    elevation = np.tile(depths, (3, 1))
    ocean_mask = np.full(elevation.shape, True)

    length_scale = _topographic_length_scale(
        elevation, ocean_mask, _DX_M, _DY_M
    )

    h1, h2 = 50.0, 10.0
    he = h1 * h2 / (h1 + h2)
    expected_km = (he * _DX_M / abs(h1 - h2)) / 1.0e3

    # both neighbors of the outlier see the same symmetric step
    assert np.isclose(length_scale[1, 2], expected_km)
    assert np.isclose(length_scale[1, 3], expected_km)
    assert np.isclose(length_scale[1, 4], expected_km)

    # far enough from the outlier, the region is locally flat
    assert np.isinf(length_scale[1, 0])


def test_length_scale_cell_width_ramps_between_thresholds():
    length_scale = np.array([0.5, 1.0, 2.5, 4.0, 10.0])
    cell_width = _length_scale_cell_width(
        length_scale=length_scale,
        max_km=3.25,
        length_scale_km=0.8125,
        threshold_low=1.0,
        threshold_high=4.0,
        exponent=1.0,
    )
    # at/below threshold_low: full refinement
    assert cell_width[0] == 0.8125
    assert cell_width[1] == 0.8125
    # at/above threshold_high: no added refinement
    assert cell_width[3] == 3.25
    assert cell_width[4] == 3.25
    # strictly between: linearly interpolated, finer than max, coarser
    # than the fully-refined value
    assert 0.8125 < cell_width[2] < 3.25
