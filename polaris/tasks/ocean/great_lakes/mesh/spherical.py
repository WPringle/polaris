import numpy as np
import xarray as xr
from geometric_features import FeatureCollection
from mpas_tools.io import open_dataset
from mpas_tools.mesh.creation.signed_distance import (
    signed_distance_from_geojson,
)
from scipy.ndimage import distance_transform_edt

from polaris.constants import get_constant
from polaris.mesh import QuasiUniformSphericalMeshStep
from polaris.mesh.spherical.coastline import signed_distance_from_ocean_mask
from polaris.mesh.spherical.graded_grid import build_graded_lon_lat


class GreatLakesSphericalMesh(QuasiUniformSphericalMeshStep):
    """
    A step for creating a JIGSAW mesh for a single Great Lake on the
    sphere, with cell width refined toward the lake's own bathymetry-
    derived shoreline

    Attributes
    ----------
    bathymetry_step : polaris.Step
        The step that prepares the lake's bathymetry

    lake_name : str
        The name of the lake
    """

    def __init__(self, component, subdir, bathymetry_step, lake_name):
        """
        Create the step

        Parameters
        ----------
        component : polaris.Component
            The component the step belongs to

        subdir : str
            The subdirectory for the step

        bathymetry_step : polaris.Step
            The step that prepares the lake's bathymetry

        lake_name : str
            The name of the lake, used to label the mesh
        """
        super().__init__(
            component=component, subdir=subdir, mesh_name=lake_name
        )
        self.bathymetry_step = bathymetry_step
        self.lake_name = lake_name

    def setup(self):
        """
        Add the bathymetry file as an input
        """
        bathy_path = self.bathymetry_step.path
        self.add_input_file(
            filename='bathymetry.nc',
            work_dir_target=f'{bathy_path}/bathymetry.nc',
        )
        super().setup()

    def build_cell_width_lat_lon(self):
        """
        Build a cell-width field that is fine (``coastal_resolution_km``)
        at the lake's bathymetry-derived shoreline, relaxing to a cap of
        ``max_resolution_km`` away from shore, and coarsening further to
        ``background_width_km`` far from the lake to keep JIGSAW
        tractable.

        Returns
        -------
        cell_width : numpy.ndarray
            m x n array of cell width in km

        lon : numpy.ndarray
            longitude in degrees (length n and between -180 and 180)

        lat : numpy.ndarray
            latitude in degrees (length m and between -90 and 90)
        """
        logger = self.logger
        section = self.config['great_lakes_mesh']

        coastal_km = section.getfloat('coastal_resolution_km')
        max_km = section.getfloat('max_resolution_km')
        transition_km = section.getfloat('coastal_transition_km')
        exponent = section.getfloat('coastal_transition_exponent')
        background_km = section.getfloat('background_width_km')
        buffer_deg = section.getfloat('region_buffer_deg')
        background_transition_km = section.getfloat('background_transition_km')
        fine_resolution_deg = section.getfloat('fine_grid_resolution_deg')
        coarse_resolution_deg = section.getfloat('coarse_grid_resolution_deg')
        growth_factor = section.getfloat('grid_growth_factor')
        slope_km = section.getfloat('slope_resolution_km')
        slope_threshold_low = section.getfloat('slope_threshold_low')
        slope_threshold_high = section.getfloat('slope_threshold_high')
        slope_exponent = section.getfloat('slope_transition_exponent')
        slope_fill_radius_km = section.getfloat('slope_fill_radius_km')
        length_scale_km = section.getfloat('length_scale_resolution_km')
        length_scale_threshold_low = section.getfloat(
            'length_scale_threshold_low_km'
        )
        length_scale_threshold_high = section.getfloat(
            'length_scale_threshold_high_km'
        )
        length_scale_exponent = section.getfloat(
            'length_scale_transition_exponent'
        )
        length_scale_fill_radius_km = section.getfloat(
            'length_scale_fill_radius_km'
        )

        logger.info(
            f'Building cell width for {self.lake_name}: refining to '
            f'{coastal_km} km at the coast, in high-slope regions, and '
            f'wherever the topographic length scale is short, plus flat '
            f'valleys enclosed by steep walls or short-length-scale '
            f'terrain (needed to resolve real glacial troughs that a '
            f'single slope threshold misses on gentler-sloped basins)'
        )

        with open_dataset('bathymetry.nc') as ds_bathy:
            bathy_lon = ds_bathy.lon.values
            bathy_lat = ds_bathy.lat.values
            elevation = ds_bathy.elevation.values

        lon_min = float(np.min(bathy_lon)) - buffer_deg
        lon_max = float(np.max(bathy_lon)) + buffer_deg
        lat_min = float(np.min(bathy_lat)) - buffer_deg
        lat_max = float(np.max(bathy_lat)) + buffer_deg

        lon, lat = build_graded_lon_lat(
            fine_lon_bounds=(lon_min, lon_max),
            fine_lat_bounds=(lat_min, lat_max),
            fine_resolution_deg=fine_resolution_deg,
            coarse_resolution_deg=coarse_resolution_deg,
            growth_factor=growth_factor,
        )

        # distance (m) to the lake's own bathymetry-derived shoreline,
        # computed on its native grid; positive over water, negative over
        # land
        ocean_mask = np.isfinite(elevation) & (elevation < 0.0)
        signed_distance = signed_distance_from_ocean_mask(
            ocean_mask, bathy_lon, bathy_lat
        )

        coastal_width_native = _coastal_cell_width(
            signed_distance=signed_distance,
            coastal_km=coastal_km,
            max_km=max_km,
            transition_km=transition_km,
            exponent=exponent,
        )

        dx_m, dy_m = _native_grid_spacing_m(bathy_lon, bathy_lat)

        slope = _bathymetric_slope(elevation, ocean_mask, dx_m, dy_m)
        slope_width_native = _slope_cell_width(
            slope=slope,
            max_km=max_km,
            slope_km=slope_km,
            threshold_low=slope_threshold_low,
            threshold_high=slope_threshold_high,
            exponent=slope_exponent,
        )

        # a trough floor fully enclosed by steep walls should be refined
        # too, even where its own local slope is too low to trigger the
        # ramp above (see the design doc's 'Slope-based refinement' section)
        closed_steep = _fill_enclosed_region(
            mask=slope >= slope_threshold_high,
            dx_m=dx_m,
            dy_m=dy_m,
            radius_km=slope_fill_radius_km,
        )
        slope_width_native = np.minimum(
            slope_width_native,
            np.where(closed_steep, slope_km, max_km),
        )

        # the topographic length scale is a scale-aware complement to the
        # raw slope threshold above (see _topographic_length_scale's
        # docstring); it is only evaluated outside the coastal-transition
        # buffer, since he collapses toward zero in shallow water near any
        # shoreline regardless of how gentle the slope actually is there
        # -- a case the coastal-distance ramp above already handles
        length_scale = _topographic_length_scale(
            elevation, ocean_mask, dx_m, dy_m
        )
        within_coastal_buffer = np.abs(signed_distance) <= (
            transition_km * 1.0e3
        )
        length_scale = np.where(within_coastal_buffer, np.inf, length_scale)

        length_scale_width_native = _length_scale_cell_width(
            length_scale=length_scale,
            max_km=max_km,
            length_scale_km=length_scale_km,
            threshold_low=length_scale_threshold_low,
            threshold_high=length_scale_threshold_high,
            exponent=length_scale_exponent,
        )

        closed_short_length_scale = _fill_enclosed_region(
            mask=length_scale <= length_scale_threshold_low,
            dx_m=dx_m,
            dy_m=dy_m,
            radius_km=length_scale_fill_radius_km,
        )
        length_scale_width_native = np.minimum(
            length_scale_width_native,
            np.where(closed_short_length_scale, length_scale_km, max_km),
        )

        coastal_width_native = np.minimum(
            coastal_width_native, slope_width_native
        )
        coastal_width_native = np.minimum(
            coastal_width_native, length_scale_width_native
        )

        cell_width = _paste_onto_graded_grid(
            coastal_width_native=coastal_width_native,
            bathy_lon=bathy_lon,
            bathy_lat=bathy_lat,
            lon=lon,
            lat=lat,
            lon_bounds=(lon_min, lon_max),
            lat_bounds=(lat_min, lat_max),
            background_km=background_km,
            fill_km=max_km,
        )

        cell_width = _coarsen_outside_region(
            cell_width=cell_width,
            lon=lon,
            lat=lat,
            lon_bounds=(lon_min, lon_max),
            lat_bounds=(lat_min, lat_max),
            background_km=background_km,
            transition_km=background_transition_km,
        )

        return cell_width, lon, lat


def _coastal_cell_width(
    signed_distance, coastal_km, max_km, transition_km, exponent
):
    """
    Blend from ``coastal_km`` at the shoreline to ``max_km`` away from it,
    symmetrically on both the water and land sides, following the same
    transition idiom used by the unified mesh's coastline refinement (see
    ``polaris.tasks.mesh.spherical.unified.sizing_field.build``).
    """
    cell_width = np.full(signed_distance.shape, max_km, dtype=float)
    transition_m = transition_km * 1.0e3
    distance = np.abs(signed_distance)

    if transition_m <= 0.0:
        at_coast = np.isclose(distance, 0.0)
        cell_width[at_coast] = coastal_km
        return cell_width

    in_transition = distance <= transition_m
    fraction = (distance[in_transition] / transition_m) ** exponent
    cell_width[in_transition] = coastal_km + fraction * (max_km - coastal_km)
    return cell_width


def _native_grid_spacing_m(bathy_lon, bathy_lat):
    """
    Native bathymetry grid spacing in meters, converting lon/lat spacing
    using the mean latitude of the basin's own bounding box (precise
    enough given a basin spans only a few degrees of latitude). Shared by
    ``_bathymetric_slope`` and ``_topographic_length_scale``, which both
    need the same native-resolution length scale.
    """
    earth_radius = get_constant('mean_radius')
    meters_per_deg = (np.pi / 180.0) * earth_radius
    dlat_deg = float(bathy_lat[1] - bathy_lat[0])
    dlon_deg = float(bathy_lon[1] - bathy_lon[0])
    mean_lat_rad = np.radians(np.mean(bathy_lat))
    dx_m = abs(dlon_deg * meters_per_deg * np.cos(mean_lat_rad))
    dy_m = abs(dlat_deg * meters_per_deg)
    return dx_m, dy_m


def _bathymetric_slope(elevation, ocean_mask, dx_m, dy_m):
    """
    Magnitude of the local bathymetric gradient (dimensionless: m of depth
    change per m of horizontal distance) on the native bathymetry grid,
    zero on land and wherever a neighboring point is land or missing (so
    it never refines the coastline itself, only interior steep features;
    the coastal-distance blend already handles the shoreline).
    """
    elevation_filled = np.where(ocean_mask, elevation, 0.0)
    dzdy, dzdx = np.gradient(elevation_filled, dy_m, dx_m)
    slope = np.hypot(dzdx, dzdy)

    # a gradient that touches a land/missing neighbor isn't a real
    # bathymetric slope; invalidate it along with points that are
    # themselves not ocean
    valid = ocean_mask.copy()
    valid[1:, :] &= ocean_mask[:-1, :]
    valid[:-1, :] &= ocean_mask[1:, :]
    valid[:, 1:] &= ocean_mask[:, :-1]
    valid[:, :-1] &= ocean_mask[:, 1:]

    return np.where(valid, slope, 0.0)


def _topographic_length_scale(elevation, ocean_mask, dx_m, dy_m):
    """
    A scale-aware alternative to raw bathymetric slope: the topographic
    length scale ``L = he / |grad(H)|``, where ``he = h1*h2/(h1+h2)`` is
    the harmonic-mean-like "effective depth" between two neighboring
    native-grid depths ``h1``, ``h2`` (the same combination that sets the
    two-layer internal wave speed ``c = sqrt(g' * he)`` in gravity-
    current/exchange-flow theory, e.g. strait sill flows). ``L`` has units
    of length: how far, at the local depth gradient, depth would need to
    change by its own effective scale -- so, unlike a raw slope threshold,
    it responds to *relative* depth change. This is what lets one shared
    threshold work across lakes with very different native steepness (see
    ``length_scale_threshold_low_km``'s comment in ``great_lakes_mesh.cfg``):
    checked directly against all four Great Lakes basins' native DEMs, one
    threshold gives similar, modest cost increases for every lake, where a
    raw slope threshold either misses Lake Michigan-Huron's real (but
    comparatively gentle) troughs or drastically over-refines Lake
    Superior.

    ``L`` is computed separately between x- and y-direction neighbor pairs
    on the native grid, each pair's value assigned to both of its
    endpoints, and the two directions combined by taking the smaller
    (more restrictive) value at each point. Returns ``inf`` (km) wherever
    two neighboring depths are equal (no local gradient) or both
    neighbors are land (``he = 0`` there), so flat or dry regions never
    drive refinement on their own.

    Returns
    -------
    length_scale_km : numpy.ndarray
        the topographic length scale in km
    """
    depth = np.where(ocean_mask, -elevation, 0.0)

    def _pair_length_scale(h1, h2, spacing_m):
        depth_sum = h1 + h2
        he = np.where(
            depth_sum > 0.0,
            h1 * h2 / np.where(depth_sum > 0.0, depth_sum, 1.0),
            0.0,
        )
        diff = np.abs(h2 - h1)
        safe_diff = np.where(diff > 0.0, diff, 1.0)
        return np.where(diff > 0.0, he * spacing_m / safe_diff, np.inf)

    length_x = _pair_length_scale(depth[:, :-1], depth[:, 1:], dx_m)
    length_y = _pair_length_scale(depth[:-1, :], depth[1:, :], dy_m)

    length_scale_x = np.full(depth.shape, np.inf)
    length_scale_x[:, :-1] = np.minimum(length_scale_x[:, :-1], length_x)
    length_scale_x[:, 1:] = np.minimum(length_scale_x[:, 1:], length_x)

    length_scale_y = np.full(depth.shape, np.inf)
    length_scale_y[:-1, :] = np.minimum(length_scale_y[:-1, :], length_y)
    length_scale_y[1:, :] = np.minimum(length_scale_y[1:, :], length_y)

    length_scale_m = np.minimum(length_scale_x, length_scale_y)
    return length_scale_m / 1.0e3


def _fill_enclosed_region(mask, dx_m, dy_m, radius_km):
    """
    Close ``mask`` via two Euclidean distance transforms -- dilate by
    thresholding the distance to the nearest ``True`` cell, then erode
    that dilated mask by thresholding the distance to its complement
    (``erosion(A) = NOT(dilation(NOT(A)))``) -- so that ``False`` pockets
    fully enclosed by a (possibly gappy) ring of ``True`` cells -- e.g. a
    flat trough floor ringed by steep trough walls, or a long-length-scale
    pocket ringed by short-length-scale walls -- are filled in too, even
    though the pocket's own cells never satisfy the ``True`` condition.
    Unlike dilating every ``True`` cell outward by ``radius_km`` alone
    (which would also expand into ordinary, unenclosed ``True`` regions
    such as a shoreline drop-off), closing only fills pockets that are
    actually enclosed, bridging wall gaps up to about ``radius_km`` wide.

    This is equivalent to morphological closing with a disk/ellipse
    structuring element of radius ``radius_km``, but computed via scipy's
    O(N) exact distance transform (which natively supports anisotropic
    ``dx_m``/``dy_m`` spacing via its ``sampling`` argument) rather than
    brute-force structuring-element morphology, whose cost scales with the
    structuring element's pixel count and so does not scale to a basin's
    full native-resolution grid.

    Shared by the slope-based and topographic-length-scale-based
    refinement criteria in ``build_cell_width_lat_lon``, which each build
    their own boolean trigger mask (``slope >= slope_threshold_high`` or
    ``length_scale <= length_scale_threshold_low_km``) and pass it here.

    Returns
    -------
    closed : numpy.ndarray of bool
        True wherever a cell is in ``mask`` itself or lies in an enclosed
        pocket that closing fills in
    """
    radius_m = radius_km * 1.0e3
    sampling = (dy_m, dx_m)

    dilated = distance_transform_edt(~mask, sampling=sampling) <= radius_m
    return distance_transform_edt(dilated, sampling=sampling) > radius_m


def _slope_cell_width(
    slope, max_km, slope_km, threshold_low, threshold_high, exponent
):
    """
    Blend from ``max_km`` to ``slope_km`` as bathymetric slope increases
    from ``threshold_low`` to ``threshold_high``, following the same
    transition idiom as ``_coastal_cell_width``.
    """
    cell_width = np.full(slope.shape, max_km, dtype=float)
    span = threshold_high - threshold_low

    if span <= 0.0:
        cell_width[slope >= threshold_high] = slope_km
        return cell_width

    in_transition = slope > threshold_low
    fraction = (
        np.clip((slope[in_transition] - threshold_low) / span, 0.0, 1.0)
        ** exponent
    )
    cell_width[in_transition] = max_km + fraction * (slope_km - max_km)
    return cell_width


def _length_scale_cell_width(
    length_scale,
    max_km,
    length_scale_km,
    threshold_low,
    threshold_high,
    exponent,
):
    """
    Blend from ``max_km`` to ``length_scale_km`` as the topographic length
    scale decreases from ``threshold_high`` to ``threshold_low`` -- the
    same transition idiom as ``_slope_cell_width``, but inverted, since a
    *smaller* length scale (not a larger one) means more refinement is
    needed.
    """
    cell_width = np.full(length_scale.shape, max_km, dtype=float)
    span = threshold_high - threshold_low

    if span <= 0.0:
        cell_width[length_scale <= threshold_low] = length_scale_km
        return cell_width

    in_transition = length_scale < threshold_high
    fraction = (
        np.clip(
            (threshold_high - length_scale[in_transition]) / span, 0.0, 1.0
        )
        ** exponent
    )
    cell_width[in_transition] = max_km + fraction * (length_scale_km - max_km)
    return cell_width


def _paste_onto_graded_grid(
    coastal_width_native,
    bathy_lon,
    bathy_lat,
    lon,
    lat,
    lon_bounds,
    lat_bounds,
    background_km,
    fill_km,
):
    """
    Interpolate the native-resolution coastal cell-width field onto the
    sub-block of the graded grid that overlaps the lake's bounding box,
    leaving ``background_km`` everywhere else.
    """
    lon_min, lon_max = lon_bounds
    lat_min, lat_max = lat_bounds

    da_native = xr.DataArray(
        coastal_width_native,
        dims=('lat', 'lon'),
        coords={'lat': bathy_lat, 'lon': bathy_lon},
    )

    fine_lon_mask = (lon >= lon_min) & (lon <= lon_max)
    fine_lat_mask = (lat >= lat_min) & (lat <= lat_max)
    fine_lon = lon[fine_lon_mask]
    fine_lat = lat[fine_lat_mask]

    # plain 1-D arrays (not sharing a dimension name) give an outer-product
    # grid interpolation, not a pointwise one
    interpolated = da_native.interp(
        lon=fine_lon, lat=fine_lat, method='linear'
    )
    interpolated = interpolated.fillna(fill_km)

    cell_width = np.full((lat.size, lon.size), background_km)
    cell_width[np.ix_(fine_lat_mask, fine_lon_mask)] = interpolated.values
    return cell_width


def _coarsen_outside_region(
    cell_width, lon, lat, lon_bounds, lat_bounds, background_km, transition_km
):
    """
    Smoothly relax ``cell_width`` toward ``background_km`` outside the
    lake's buffered bounding box, following the same signed-distance blend
    used by ``isomip_plus``'s small-region spherical mesh (see
    ``polaris.tasks.ocean.isomip_plus.mesh.spherical``), which confirms
    ``signed_distance_from_geojson`` works on the full sphere, not just a
    plane.
    """
    earth_radius = get_constant('mean_radius')
    fc = _bounding_box_feature(lon_bounds, lat_bounds)
    region_signed_distance = signed_distance_from_geojson(
        fc, lon, lat, earth_radius, max_length=0.25
    )

    trans_width_m = transition_km * 1.0e3
    weights = np.clip(region_signed_distance / trans_width_m, 0.0, 1.0)
    return cell_width * (1.0 - weights) + background_km * weights


def _bounding_box_feature(lon_bounds, lat_bounds):
    """
    Build a ``geometric_features.FeatureCollection`` polygon for a lon/lat
    bounding box, for use with ``signed_distance_from_geojson``.
    """
    lon_min, lon_max = lon_bounds
    lat_min, lat_max = lat_bounds
    coordinates = [
        [
            [lon_min, lat_min],
            [lon_max, lat_min],
            [lon_max, lat_max],
            [lon_min, lat_max],
            [lon_min, lat_min],
        ]
    ]
    features = [
        {
            'type': 'Feature',
            'properties': {
                'name': 'Great Lakes mesh region',
                'component': 'ocean',
                'object': 'region',
                'author': 'Polaris',
            },
            'geometry': {'type': 'Polygon', 'coordinates': coordinates},
        }
    ]
    return FeatureCollection(features=features)
