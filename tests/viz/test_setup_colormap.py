"""
Unit tests for :py:func:`polaris.viz.spherical.setup_colormap`.
"""

from polaris.config import PolarisConfigParser
from polaris.viz.spherical import setup_colormap


def _config(colormap_levels=None):
    config = PolarisConfigParser()
    config.add_section('viz')
    config.set('viz', 'colormap_name', 'viridis')
    config.set('viz', 'norm_type', 'linear')
    config.set('viz', 'norm_args', '{"vmin": 0.0, "vmax": 1.0}')
    if colormap_levels is not None:
        config.set('viz', 'colormap_levels', str(colormap_levels))
    return config


def test_colormap_levels_resamples_the_colormap():
    colormap, _, _ = setup_colormap(_config(colormap_levels=12), 'viz')
    assert colormap.N == 12


def test_without_colormap_levels_the_colormap_is_not_resampled():
    colormap, _, _ = setup_colormap(_config(), 'viz')
    assert colormap.N == 256


def test_colormap_levels_keeps_under_and_over_colors():
    config = _config(colormap_levels=12)
    config.set('viz', 'under_color', 'black')
    config.set('viz', 'over_color', 'white')

    colormap, _, _ = setup_colormap(config, 'viz')

    assert colormap.N == 12
    assert tuple(colormap.get_under()) == (0.0, 0.0, 0.0, 1.0)
    assert tuple(colormap.get_over()) == (1.0, 1.0, 1.0, 1.0)
