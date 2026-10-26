import numpy as np

from polaris.mesh.spherical.graded_grid import build_graded_lon_lat


def test_graded_grid_is_monotonic_and_spans_the_globe():
    lon, lat = build_graded_lon_lat(
        fine_lon_bounds=(-93.0, -75.0),
        fine_lat_bounds=(40.0, 50.0),
        fine_resolution_deg=0.5,
        coarse_resolution_deg=5.0,
        growth_factor=1.3,
    )
    assert np.all(np.diff(lon) > 0.0)
    assert np.all(np.diff(lat) > 0.0)
    assert np.isclose(lon[0], -180.0)
    assert np.isclose(lon[-1], 180.0)
    assert np.isclose(lat[0], -90.0)
    assert np.isclose(lat[-1], 90.0)
    assert np.diff(lon).max() <= 5.0 + 1.0e-10
    assert np.diff(lat).max() <= 5.0 + 1.0e-10


def test_graded_grid_is_uniform_and_fine_within_the_fine_region():
    lon, _ = build_graded_lon_lat(
        fine_lon_bounds=(-93.0, -75.0),
        fine_lat_bounds=(40.0, 50.0),
        fine_resolution_deg=0.5,
        coarse_resolution_deg=5.0,
        growth_factor=1.3,
    )
    in_fine_region = (lon >= -93.0) & (lon <= -75.0)
    fine_diffs = np.diff(lon[in_fine_region])
    assert np.allclose(fine_diffs, 0.5)


def test_graded_grid_steps_grow_geometrically_outside_the_fine_region():
    lon, _ = build_graded_lon_lat(
        fine_lon_bounds=(-93.0, -75.0),
        fine_lat_bounds=(40.0, 50.0),
        fine_resolution_deg=0.5,
        coarse_resolution_deg=5.0,
        growth_factor=1.3,
    )
    diffs = np.diff(lon)
    idx_fine_end = np.searchsorted(lon, -75.0)
    # the step sizes leaving the fine region should grow monotonically
    # until they saturate at the coarse resolution
    growing = diffs[idx_fine_end : idx_fine_end + 8]
    assert np.all(np.diff(growing) > 0.0)
    assert np.isclose(growing[0], 0.5 * 1.3)


def test_graded_grid_reduces_to_uniform_when_resolutions_match():
    lon, lat = build_graded_lon_lat(
        fine_lon_bounds=(-10.0, 10.0),
        fine_lat_bounds=(-5.0, 5.0),
        fine_resolution_deg=1.0,
        coarse_resolution_deg=1.0,
        growth_factor=1.5,
    )
    assert np.allclose(np.diff(lon), 1.0)
    assert np.allclose(np.diff(lat), 1.0)
