import numpy as np


def build_graded_lon_lat(
    fine_lon_bounds,
    fine_lat_bounds,
    fine_resolution_deg,
    coarse_resolution_deg,
    growth_factor=1.15,
):
    """
    Build a monotonic, non-uniform longitude and latitude grid that is fine
    within a region of interest and coarsens geometrically outward to a
    coarse background resolution.

    ``QuasiUniformSphericalMeshStep.make_jigsaw_mesh()`` hands JIGSAW's
    ``ELLIPSOID-GRID`` sizing function two 1-D coordinate arrays that only
    need to be monotonic, not uniformly spaced.  A uniform global grid fine
    enough to resolve a small regional coastline would have tens to hundreds
    of millions of points; grading the grid so only the region of interest
    is fine keeps the point count tractable.

    Parameters
    ----------
    fine_lon_bounds : tuple of float
        The (min, max) longitude in degrees of the fine region, in -180 to
        180

    fine_lat_bounds : tuple of float
        The (min, max) latitude in degrees of the fine region, in -90 to 90

    fine_resolution_deg : float
        The grid spacing in degrees within the fine region

    coarse_resolution_deg : float
        The grid spacing in degrees far from the fine region

    growth_factor : float, optional
        The factor by which grid spacing grows at each step outward from
        the fine region until it reaches ``coarse_resolution_deg``

    Returns
    -------
    lon : numpy.ndarray
        A monotonically increasing 1-D longitude array spanning -180 to 180

    lat : numpy.ndarray
        A monotonically increasing 1-D latitude array spanning -90 to 90
    """
    lon = _graded_axis(
        fine_min=fine_lon_bounds[0],
        fine_max=fine_lon_bounds[1],
        domain_min=-180.0,
        domain_max=180.0,
        fine_resolution=fine_resolution_deg,
        coarse_resolution=coarse_resolution_deg,
        growth_factor=growth_factor,
    )
    lat = _graded_axis(
        fine_min=fine_lat_bounds[0],
        fine_max=fine_lat_bounds[1],
        domain_min=-90.0,
        domain_max=90.0,
        fine_resolution=fine_resolution_deg,
        coarse_resolution=coarse_resolution_deg,
        growth_factor=growth_factor,
    )
    return lon, lat


def _graded_axis(
    fine_min,
    fine_max,
    domain_min,
    domain_max,
    fine_resolution,
    coarse_resolution,
    growth_factor,
):
    """
    Build one monotonic, graded coordinate axis: uniform at
    ``fine_resolution`` between ``fine_min`` and ``fine_max``, growing
    geometrically out to ``coarse_resolution`` in each direction, then
    uniform at ``coarse_resolution`` out to the domain bounds.
    """
    if fine_min < domain_min or fine_max > domain_max:
        raise ValueError(
            'Fine bounds must lie within the domain bounds: '
            f'({fine_min}, {fine_max}) not in ({domain_min}, {domain_max})'
        )
    if fine_min >= fine_max:
        raise ValueError('fine_min must be less than fine_max')
    if fine_resolution <= 0.0 or coarse_resolution <= 0.0:
        raise ValueError('Resolutions must be positive')
    if coarse_resolution < fine_resolution:
        raise ValueError(
            'coarse_resolution must be greater than or equal to '
            'fine_resolution'
        )
    if growth_factor <= 1.0:
        raise ValueError('growth_factor must be greater than one')

    fine_count = int(np.ceil((fine_max - fine_min) / fine_resolution)) + 1
    fine_points = np.linspace(fine_min, fine_max, fine_count)

    left_points = _grow_points(
        start=fine_points[0],
        end=domain_min,
        first_step=fine_resolution,
        coarse_step=coarse_resolution,
        growth_factor=growth_factor,
        direction=-1.0,
    )
    right_points = _grow_points(
        start=fine_points[-1],
        end=domain_max,
        first_step=fine_resolution,
        coarse_step=coarse_resolution,
        growth_factor=growth_factor,
        direction=1.0,
    )

    axis = np.concatenate([left_points[::-1], fine_points, right_points])
    return axis


def _grow_points(
    start, end, first_step, coarse_step, growth_factor, direction
):
    """
    Walk outward from ``start`` toward ``end`` (``direction`` is +1 or -1),
    growing the step size geometrically from ``first_step`` up to
    ``coarse_step``, then continuing at ``coarse_step``.  Returns the new
    points in the order they were generated (increasing distance from
    ``start``), not including ``start`` itself.  The final point is always
    exactly ``end``.
    """
    points = []
    position = start
    step = first_step
    remaining = direction * (end - start)
    while remaining > 0.0:
        step = min(step * growth_factor, coarse_step)
        if step >= remaining:
            position = end
            points.append(position)
            break
        position += direction * step
        points.append(position)
        remaining -= step

    return np.array(points)
