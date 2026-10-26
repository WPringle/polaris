import numpy as np

from polaris.tasks.ocean.great_lakes.mesh.spherical import (
    _fill_enclosed_steep_regions,
)

_DX_M = 100.0
_DY_M = 100.0
_THRESHOLD_HIGH = 0.1
_BOX = (16, 25, 16, 25)
_CENTER = (20, 20)


def _box_ring_slope(shape, box, gap=None, gap_width=1):
    """
    A slope array that is ``_THRESHOLD_HIGH`` (steep) along the perimeter
    of ``box`` (a (row_start, row_stop, col_start, col_stop) tuple) and
    zero everywhere else, with ``gap_width`` cells of the top edge cleared
    starting at column ``gap`` if given, to test that closing can bridge a
    small gap in an otherwise-enclosing wall.
    """
    slope = np.zeros(shape)
    row_start, row_stop, col_start, col_stop = box
    slope[row_start, col_start:col_stop] = _THRESHOLD_HIGH
    slope[row_stop - 1, col_start:col_stop] = _THRESHOLD_HIGH
    slope[row_start:row_stop, col_start] = _THRESHOLD_HIGH
    slope[row_start:row_stop, col_stop - 1] = _THRESHOLD_HIGH
    if gap is not None:
        slope[row_start, gap : gap + gap_width] = 0.0
    return slope


def test_fill_enclosed_steep_regions_fills_a_fully_enclosed_pocket():
    slope = _box_ring_slope(shape=(41, 41), box=_BOX)

    # too small a radius leaves the interior pocket unfilled
    closed_small = _fill_enclosed_steep_regions(
        slope=slope,
        threshold_high=_THRESHOLD_HIGH,
        dx_m=_DX_M,
        dy_m=_DY_M,
        radius_km=0.2,
    )
    assert not closed_small[_CENTER]

    # a large-enough radius fills the low-slope interior, fully ringed by
    # steep cells, without spilling out far beyond the ring
    closed_large = _fill_enclosed_steep_regions(
        slope=slope,
        threshold_high=_THRESHOLD_HIGH,
        dx_m=_DX_M,
        dy_m=_DY_M,
        radius_km=0.4,
    )
    assert closed_large[_CENTER]
    assert not closed_large[0, 0]
    assert not closed_large[40, 40]


def test_fill_enclosed_steep_regions_bridges_a_gap_in_the_wall():
    slope = _box_ring_slope(shape=(41, 41), box=_BOX, gap=19, gap_width=3)
    closed = _fill_enclosed_steep_regions(
        slope=slope,
        threshold_high=_THRESHOLD_HIGH,
        dx_m=_DX_M,
        dy_m=_DY_M,
        radius_km=0.4,
    )
    assert closed[_CENTER]


def test_fill_enclosed_steep_regions_does_not_fill_an_open_steep_line():
    # a straight, open line of steep cells (e.g. an ordinary shoreline
    # drop-off) encloses nothing, so closing should not balloon out to
    # fill either side of it
    slope = np.zeros((41, 41))
    slope[20, :] = _THRESHOLD_HIGH
    closed = _fill_enclosed_steep_regions(
        slope=slope,
        threshold_high=_THRESHOLD_HIGH,
        dx_m=_DX_M,
        dy_m=_DY_M,
        radius_km=0.4,
    )
    assert not closed[0, 20]
    assert not closed[40, 20]
