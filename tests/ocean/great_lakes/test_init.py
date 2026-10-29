"""
Unit tests for :py:mod:`polaris.tasks.ocean.great_lakes.init`.

Self-contained: no file I/O and no full Polaris step framework, testing
only the pure functions the step's ``run()`` calls.
"""

import numpy as np
import xarray as xr

from polaris.tasks.ocean.great_lakes.init import (
    _bottom_depth_from_elevation,
    _idealized_temperature,
    _salinity_by_basin,
)

# ---------------------------------------------------------------------------
# _bottom_depth_from_elevation
# ---------------------------------------------------------------------------


def test_bottom_depth_from_elevation_converts_sign():
    elevation = xr.DataArray([-10.0, -100.0, -5.0])
    depth = _bottom_depth_from_elevation(elevation, min_bottom_depth=1.0)
    assert np.allclose(depth.values, [10.0, 100.0, 5.0])


def test_bottom_depth_from_elevation_clamps_shallow_cells():
    elevation = xr.DataArray([-0.5, -2.0, -10.0])
    depth = _bottom_depth_from_elevation(elevation, min_bottom_depth=3.0)
    assert np.allclose(depth.values, [3.0, 3.0, 10.0])


def test_bottom_depth_from_elevation_clamps_fill_value_cells():
    elevation = xr.DataArray([-10.0, np.nan, -20.0])
    depth = _bottom_depth_from_elevation(elevation, min_bottom_depth=5.0)
    assert np.allclose(depth.values, [10.0, 5.0, 20.0])
    assert not np.any(np.isnan(depth.values))


# ---------------------------------------------------------------------------
# _idealized_temperature
# ---------------------------------------------------------------------------


def test_idealized_temperature_matches_surface_and_bottom_asymptotes():
    # far above and far below the thermocline, temperature should sit at
    # (within floating-point reach of) the two end-member values
    z_mid = xr.DataArray([[-0.1, -500.0]], dims=('nCells', 'nVertLevels'))
    temperature = _idealized_temperature(
        z_mid,
        surface_temperature=20.0,
        bottom_temperature=4.0,
        thermocline_depth=15.0,
        thermocline_width=5.0,
    )
    assert temperature.dims == ('Time', 'nCells', 'nVertLevels')
    assert np.isclose(temperature.values[0, 0, 0], 20.0, atol=0.1)
    assert np.isclose(temperature.values[0, 0, 1], 4.0, atol=1e-6)


def test_idealized_temperature_is_monotonic_with_depth():
    z_mid = xr.DataArray(
        [np.linspace(-0.1, -200.0, 50)], dims=('nCells', 'nVertLevels')
    )
    temperature = _idealized_temperature(
        z_mid,
        surface_temperature=20.0,
        bottom_temperature=4.0,
        thermocline_depth=15.0,
        thermocline_width=5.0,
    )
    values = temperature.values[0, 0, :]
    assert np.all(np.diff(values) <= 0.0)


# ---------------------------------------------------------------------------
# _salinity_by_basin
# ---------------------------------------------------------------------------


def test_salinity_by_basin_varies_by_basin_not_by_level():
    temperature = xr.DataArray(
        np.zeros((1, 3, 4)), dims=('Time', 'nCells', 'nVertLevels')
    )
    basin_index = xr.DataArray([0, 1, 0], dims=('nCells',))
    salinity = _salinity_by_basin(
        temperature, basin_index, salinity_by_basin=[0.05, 0.2]
    )

    assert salinity.dims == temperature.dims
    # basin 0's two cells both get 0.05, basin 1's cell gets 0.2, uniform
    # across all vertical levels
    assert np.allclose(salinity.values[0, 0, :], 0.05)
    assert np.allclose(salinity.values[0, 1, :], 0.2)
    assert np.allclose(salinity.values[0, 2, :], 0.05)
