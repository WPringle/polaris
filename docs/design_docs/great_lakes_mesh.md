# Great Lakes regional mesh

Creation date: 2026/10/22

Contributors: Claude Code (with wpringle)

## Summary

This adds a Polaris task that builds an Omega/MPAS-Ocean mesh for a single
Laurentian Great Lake from a bathymetry DEM, with resolution refined toward
the lake's own shoreline. No existing Polaris mesh pipeline fits: the
unified-mesh JIGSAW pipeline (`polaris/tasks/mesh/spherical/unified/`)
builds a single global mesh from shared global coastline and river
datasets, and cannot be scoped to a small region. `isomip_plus`'s
small-region spherical mesh uses a synthetic coastline, not a real DEM.

Success is a culled MPAS mesh per basin with cell width 812.5 m at the
shoreline, capped at 3.25 km in the open lake, and `base_elevation` from the
DEM on the mesh.

Each basin is built and culled as its own closed sub-domain. Lake Michigan
and Lake Huron are combined into one basin (`michigan_huron`), since they
are hydrologically a single lake connected by the Straits of Mackinac and
their downloaded DEMs tile together exactly; the remaining basins
(Superior, Michigan-Huron, Erie, Ontario) are not connected to each other
or to the St. Lawrence River in this phase; see [Decisions](#decisions).

## Requirements

### Requirement: Mesh resolution

Date last modified: 2026/10/26

Cell width is 0.8125 km at the shoreline, relaxing to a cap of 3.25 km away
from shore, on both the water and land sides of the coastline. 0.8125 km
is `max_resolution_km / 4` and close to 3x the ~270 m pixel spacing of the
coarsened bathymetry used for the final remap, so the remap source stays
several times finer than the target resolution.

Cell width is also 0.8125 km (812.5 m) -- matching the shoreline target --
wherever bathymetric slope exceeds `slope_threshold_high` = 0.1
(dimensionless: m of depth change per m of horizontal distance, i.e. a 10%
grade), ramping in from `max_resolution_km` starting at
`slope_threshold_low` = 0.02 (a 2% grade); see "Algorithm Design:
Slope-based refinement" for why and how this is computed. The same
0.8125 km resolution is also reached wherever a second, scale-aware
criterion -- a topographic length scale `L = he / |grad(H)|` -- drops to
`length_scale_threshold_low_km` = 1.5 km or below, ramping in from
`length_scale_threshold_high_km` = 4.0 km; see "Algorithm Design:
Topographic length-scale refinement" for why this second criterion is
needed alongside the slope-based one.

### Requirement: Mesh domain

Date last modified: 2026/10/23

The mesh covers one basin's wet cells, as determined by its own bathymetry
(which may be merged from more than one downloaded DEM, for
`michigan_huron`). Cells far from the basin are culled regardless of
depth.

### Requirement: Bathymetry on the mesh

Date last modified: 2026/10/22

The culled mesh carries `base_elevation`, conservatively remapped from the
source DEM, consistent with the `base_elevation` naming already used by
`polaris/tasks/e3sm/init/topo/combine`. For a cell with less than 100%
ocean fraction (survived culling at `min_ocean_fraction` but still
partly land), `base_elevation` is the area-weighted average depth of the
cell's wet portion only -- not diluted by land elevation. See "Algorithm
Design: Topography normalization" below.

## Algorithm Design

### Algorithm Design: Cell-width field

Date last modified: 2026/10/22

`QuasiUniformSphericalMeshStep.make_jigsaw_mesh()` writes JIGSAW's sizing
function as an `ELLIPSOID-GRID` hfun from two 1-D lon/lat coordinate
arrays, which must be monotonic but need not be uniform. A uniform global
grid fine enough for 812.5 m resolution would have on the order of 10^8
points. Instead, `polaris.mesh.spherical.graded_grid.build_graded_lon_lat`
builds a grid that is uniform and fine over the lake's bounding box
(padded by `region_buffer_deg`) and coarsens geometrically outward to a
background resolution.

Within the fine region, cell width is computed on the DEM's own native
grid: `polaris.mesh.spherical.coastline.signed_distance_from_ocean_mask`
gives signed distance to the DEM-derived shoreline (`elevation < 0`), and
the same blend used by the unified mesh's coastline refinement
(`_apply_transition` in
`polaris/tasks/mesh/spherical/unified/sizing_field/build.py`) ramps cell
width from `coastal_resolution_km` at distance zero to `max_resolution_km`
at `coastal_transition_km`, on the unsigned distance (so land cells near
the shore are refined too, for a clean post-cull coastline). A second,
independent refinement based on bathymetric slope (see "Algorithm Design:
Slope-based refinement" below) is combined with this coastal field by
taking the finer (smaller) of the two at each native grid point. The
combined field is then interpolated onto the graded grid's fine
sub-block.

### Algorithm Design: Slope-based refinement

Date last modified: 2026/10/26

Coastal-distance refinement alone leaves every offshore cell at
`max_resolution_km` regardless of what's underneath it, so a steep,
narrow trench far from shore gets the same coarse cells as flat open
water, and conservative remapping onto cells that wide averages the
trench's extreme depths away. Checking Lake Superior's culled mesh after
its first end-to-end run showed exactly this: `base_elevation` only
reached -350 m against a native-DEM deepest point of -405 m (the lake's
known actual maximum depth is close to the latter).

`_bathymetric_slope` computes `|grad(elevation)|` (dimensionless, m of
depth change per m of horizontal distance) on the native bathymetry grid
via `np.gradient`, converting lon/lat spacing to meters with the mean
latitude of the basin's own bounding box (precise enough given a basin
spans only a few degrees of latitude). It is defined only where a point
and all four of its grid neighbors are ocean, so land and
missing-bathymetry edges never trigger slope-based refinement -- the
coastal-distance blend already owns the shoreline. `_slope_cell_width`
then blends from `max_resolution_km` to `slope_resolution_km` as slope
rises from `slope_threshold_low` to `slope_threshold_high`, following the
same transition idiom as the coastal blend, and the result is combined
with the coastal-distance field via `np.minimum` (finest wins).

Threshold defaults (`slope_threshold_low = 0.02`,
`slope_threshold_high = 0.1`) come from checking the actual slope
distribution of Lake Superior's and Lake Erie's native DEMs: median
water-cell slope is around 0.006-0.03 and the 90th percentile around
0.03-0.06, so 0.1 targets roughly the steepest 1% of the bathymetry (the
slope near Superior's deepest point reaches ~0.23), not ordinary
basin-floor roughness. With these defaults, about 18% of Superior's water
cells get some refinement and about 1% are refined all the way to
`slope_resolution_km` (defaults to `coastal_resolution_km`).

With slope refinement alone, Superior's end-to-end run still only reached
-367.5 m at its deepest cell (46.91N, 86.60W), over the same native-DEM
pixel at -405.1 m used above. The reason: `_bathymetric_slope` /
`_slope_cell_width` act on each cell's own local slope, and the slope
*at the deepest pixel itself* is only about 0.004, well under
`slope_threshold_low`. The steep walls surrounding the trough (slope up
to 0.23) do get refined, but a flat trough floor ringed by steep walls
never qualifies on its own slope, no matter how steep its surroundings
are. Simulating conservative-remap area-averaging at this exact trough,
using the same 3x-coarsened source the remap step reads, confirmed
resolution would otherwise fix this: a 0.925 km cell there reaches about
-399 m (vs. -355 m at 3.25 km, and -370 m at the 2.83 km the trough
actually got) -- so the fix only needs to make `slope_resolution_km`
cell width actually reach the trough floor, not add a new, finer
resolution tier.

`_fill_enclosed_steep_regions` closes this gap: it morphologically closes
the mask of cells with `slope >= slope_threshold_high`, so low-slope
pockets fully enclosed by a ring of steep cells -- such as a trough floor
-- are filled in too, without expanding refinement along *open*,
unenclosed steep terrain such as an ordinary shoreline drop-off, which a
simple fixed-radius dilation from every steep cell would do. The closed
mask is then forced to `slope_resolution_km` and combined into the
slope-width field via `np.minimum`, same as the ramp above.

The closing is computed via two Euclidean distance transforms
(`scipy.ndimage.distance_transform_edt`, with `sampling` set to the native
grid's anisotropic pixel spacing) rather than
`scipy.ndimage.binary_closing` with an explicit structuring element: dilate
by thresholding the distance to the nearest steep cell, then erode that
dilated mask by thresholding the distance to its complement, using the
identity `erosion(A) = NOT(dilation(NOT(A)))`. An explicit-structuring-
element closing with a ~2 km-radius ellipse did not finish within two
minutes on Superior's real ~41M-point native grid (4201 x 9841), since its
cost scales with the structuring element's pixel count; the two-EDT
version is mathematically equivalent (verified to reproduce identical
results on both synthetic enclosed-pocket/gapped-wall/open-line cases and
the real data) but completes in about 8.8 s, since each
`distance_transform_edt` call is O(N) regardless of radius.

A plain radius-based dilation was ruled out on cost: even a 0.5 km
dilation radius from every `slope >= slope_threshold_high` cell would add
about 2,038 km^2 lake-wide (on top of the ~755 km^2 already steep),
because many ordinary shoreline/shoal areas are steep too, not just
troughs. Calibrating against Superior's southeast trough instead: its
steep "wall" has small gaps, and closing needs a radius of about 1.85 km
to bridge them and fully enclose the floor. At `slope_fill_radius_km`'s
default of 2.0 km (just above that calibrated minimum), the lake-wide area
newly refined is about 1,319 km^2 (~1.6% of the lake) -- still a small
fraction of what blind dilation would cost at the same radius
(~8,600 km^2).

Outside the buffered bounding box, cell width relaxes further to
`background_width_km` over `background_transition_km`, using
`signed_distance_from_geojson` against the box polygon -- the same
primitive `isomip_plus`'s `SphericalMesh` uses, confirming it works on the
full sphere, not only a plane.

### Algorithm Design: Topographic length-scale refinement

Date last modified: 2026/10/28

The "Testing and Validation: End-to-end mesh" section's original claim
that "only Superior's bathymetry has an enclosed steep-walled trough" was
based on checking only each lake's single globally-deepest cell, not
scanning the whole basin. Recomputing `_bathymetric_slope` directly
against Michigan-Huron's native `bathymetry.nc` (same function, same
production code path) showed two things that claim missed:

- Michigan-Huron's native DEM is a genuinely different population from
  Superior's/Erie's, not just a lake with a flatter deepest point: its
  water-cell slope is roughly half of Superior's at *every* percentile
  checked (median 0.0039 vs. 0.006, 90th 0.017 vs. 0.034, 99th 0.053 vs.
  0.097). `slope_threshold_high = 0.1` was calibrated only against
  Superior's and Erie's distributions (see "Algorithm Design: Slope-based
  refinement"), so applied unchanged to Michigan-Huron it sits above
  nearly all of its real bathymetric relief, not just ordinary
  basin-floor roughness. Connected-component labeling of water cells
  below -100 m finds several glacially-scoured Georgian Bay pockets (true
  depth 108-126 m) whose own bounding-wall slope tops out at 0.04-0.06 --
  real relief, but below `slope_threshold_high` -- so
  `_fill_enclosed_steep_regions` never closes around them and they sit at
  `max_resolution_km` almost throughout. Checked against the actual
  end-to-end Michigan-Huron run
  (`topo_remap_culled/topography_remapped.nc`): a pocket near
  (-80.6, 44.9), true minimum -113.5 m, was smoothed to only -89.2 m (a
  24.3 m / 21% bias), and a pocket near (-80.8, 44.9), true minimum
  -123.4 m, reached only -118.0 m (5.4 m bias).
- A separate, undocumented instance of the same class of bug turned up in
  Lake Ontario, for a different reason: its deepest point (-243.9 m true)
  was recorded as only -225.7 m (18.2 m bias), not because its wall slope
  was too gentle, but because `slope_fill_radius_km = 2.0` falls just
  short of the radius needed to bridge a gap in its trough wall.

Lowering `slope_threshold_low`/`slope_threshold_high` and raising
`slope_fill_radius_km` for Michigan-Huron alone (checked directly against
an end-to-end rerun with `slope_threshold_low = 0.01`,
`slope_threshold_high = 0.04`, `slope_fill_radius_km = 3.0`) did fix the
Georgian Bay pockets (33,884 -> 46,042 cells, the -113.5 m pocket's bias
down from 24.3 m to 2.6 m), but a single uniform *absolute* slope
threshold cannot serve all four lakes at once: lowered enough to catch
Michigan-Huron's gentle troughs, the same threshold applied to Superior
pushes its closed-and-refined fraction from 2.3% to about 17% of the lake
(+86% cells, 23,433 -> ~43,500) for no corresponding benefit, since
Superior's own trough is already well-captured at the existing threshold.
A per-basin override of the slope options would fix Michigan-Huron and
Ontario without over-refining Superior/Erie, but at the cost of the four
lakes no longer sharing one config -- rejected in favor of a refinement
criterion that responds to *relative* depth change instead of an
absolute slope value, so one shared threshold can work across lakes with
very different native steepness.

`_topographic_length_scale` computes, at each native-grid point, a
topographic length scale `L = he / |grad(H)|`, where
`he = h1*h2/(h1+h2)` is the harmonic-mean-like "effective depth" between
two depths `h1`, `h2` on either side of a grid step. `he` is the same
combination that sets the two-layer internal wave speed
`c = sqrt(g' * he)` in gravity-current/exchange-flow theory (e.g. strait
sill flows): a large, gentle basin and a steep-but-shallow shelf can have
the same raw slope but very different `he`, and so very different `L`.
`L` has units of length -- how far, at the local depth gradient, depth
would need to change by its own effective scale -- and is computed
pairwise along each grid axis (`x` and `y` separately, using the same
lon/lat-to-meters conversion `_bathymetric_slope` already used, now
factored into a shared `_native_grid_spacing_m` helper used by both),
then combined via an elementwise minimum, so each point's `L` reflects
its steepest nearby relative change in either direction.
`_length_scale_cell_width` then ramps resolution from `max_resolution_km`
at `length_scale_threshold_high_km` down to `length_scale_resolution_km`
at `length_scale_threshold_low_km` -- the mirror image of
`_slope_cell_width`'s ramp, since a *small* `L` means more refinement.

`he` collapses toward zero in shallow water next to any shoreline,
regardless of how gentle the real underwater slope is there, which would
otherwise make `L` spuriously small along most of the coastline -- a
near-coast artifact confirmed directly against both synthetic and real
bathymetry (`L` dropping to 0 immediately next to land even where the
true offshore slope is flat; see
`test_topographic_length_scale_is_zero_next_to_land`). Since the
coastal-distance ramp (see "Algorithm Design: Cell-width field") already
refines this entire near-shore band on its own, `build_cell_width_lat_lon`
sets `L = inf` wherever a point falls within `coastal_transition_km` of
the shoreline, reusing the exact `signed_distance` array already computed
once per lake via `signed_distance_from_ocean_mask` for that ramp (no new
distance computation needed). Flat, long-`L` trough floors fully enclosed
by a ring of short-`L` walls are closed the same way steep-walled trough
floors are for the slope path: `_fill_enclosed_steep_regions` was
generalized into `_fill_enclosed_region(mask, dx_m, dy_m, radius_km)`,
taking a precomputed boolean trigger mask (`slope >= slope_threshold_high`
for the existing path; `L <= length_scale_threshold_low_km` outside the
coastal buffer for the new one) instead of computing it internally, so
both paths share the same two-EDT morphological-closing implementation.

Checked directly against all four lakes' native DEMs, with the
coastal-buffer exclusion in place, one shared parameter set
(`length_scale_threshold_low_km = 1.5`, `length_scale_threshold_high_km =
4.0`, `length_scale_fill_radius_km = 3.0`) gives similar, modest cost
increases for every lake and fixes both known bugs without per-lake
tuning. End-to-end reruns of all four basins (Superior and Michigan-Huron
via a `shared`-QOS Perlmutter batch job, account `m3780`; Erie and
Ontario interactively; all with `MPICH_GPU_SUPPORT_ENABLED=0`) confirm
the diagnostic -- see "Testing and Validation: End-to-end mesh" for the
full results table.

Unit tests for the new functions
(`tests/ocean/great_lakes/mesh/test_spherical.py`:
`test_topographic_length_scale_flat_region_is_infinite`,
`test_topographic_length_scale_is_zero_next_to_land`,
`test_topographic_length_scale_matches_pairwise_formula`,
`test_length_scale_cell_width_ramps_between_thresholds`, plus three
`_fill_enclosed_region` cases shared with the slope path) exercise the
pairwise formula, the near-coast artifact, the ramp, and the
closing/bridging/non-filling behavior on synthetic grids.

### Algorithm Design: Merging Michigan and Huron's bathymetry

Date last modified: 2026/10/23

Lake Michigan and Lake Huron share the same water-surface elevation and
are connected by the ~5 km wide Straits of Mackinac, so they are
hydrologically one lake, even though NOAA ships their bathymetry as two
separate per-lake DEMs. Checking the two downloaded files directly showed
they tile together exactly: identical 3 arc-second grid spacing, and the
shared boundary column at `lon = -84.5` is bit-for-bit identical between
the two files (NOAA evidently split one larger survey into two per-lake
products with a one-column overlap). A transect through the strait at
`lat = 45.8` is continuously wet across the join with zero gaps.

`PrepareBathymetryStep` therefore accepts a list of source files
(`bathymetry_filenames`, comma- or whitespace-separated) rather than a
single path. When a basin has more than one source,
`_merge_bathymetry_grids` places each source's elevation array into its
own slice of an output grid sized to the union of all sources' lon/lat
extents (verifying all sources share the same grid spacing first), with
NaN outside the union of their footprints -- cell classification already
gates on `np.isfinite(elevation)` before testing `elevation < 0`
(`GreatLakesSphericalMesh.build_cell_width_lat_lon`,
`WetMaskOnMeshStep.run`), so merged-grid NaN is handled the same as land
with no further changes needed downstream. Every step after
`PrepareBathymetryStep` (mesh building, culling, remapping) is already
basin-agnostic -- none of them assume a single source file or a
particular lake's geometry -- so combining a basin from more than one DEM
required no changes outside `bathymetry.py` and the task/lake list.

### Algorithm Design: Culling

Date last modified: 2026/10/22

`isomip_plus`'s `TopoMap`/`TopoRemap`/`CullMesh` sequence runs a
conservative remap of the full DEM onto the mesh twice: once pre-cull to
get a wet/dry mask for culling, once post-cull for the final
`base_elevation`/`oceanFraction` fields. For this task the pre-cull mesh
is the uncull base mesh (tens of thousands of cells, most of which are
discarded by culling), so a full conservative remap against the DEM's
native grid there is wasted work just to decide which cells to keep.

`WetMaskOnMeshStep` replaces that first pass with cheap nearest-neighbor
point sampling of the DEM at each base-mesh cell's center (gated to the
DEM's own lon/lat extent, since `.sel(..., method='nearest')` otherwise
nearest-matches cells far outside the lake to an edge pixel), followed by
`polaris.mesh.connectivity.connected_to_seeds`/`seed_mask_from_points` to
drop isolated single-cell "islands" that a coarse graded-grid cell's
center can spuriously land on near the coastal-transition zone. The
result seeds `GreatLakesCullMesh`, which culls cells below
`min_ocean_fraction` via `mpas_tools.mesh.conversion.cull()`, reusing
`isomip_plus`'s `_land_mask_from_topo` idiom directly.

The full conservative remap (`GreatLakesTopoMap`/`GreatLakesTopoRemap`)
then runs only once, on the already-culled mesh, for the final
`base_elevation`/`oceanFraction` fields -- where an accurate
area-weighted ocean fraction actually matters.

This ordering (topography before cull) differs from
`polaris/tasks/e3sm/init/topo/cull`, where land/ocean masking is built
independently of mesh resolution; here resolution itself depends on the
DEM's shoreline, so the DEM must be processed first regardless.

### Algorithm Design: Topography normalization

Date last modified: 2026/10/23

A conservative remap of raw elevation (positive on land, negative in
water) onto a cell with wet fraction `f` gives
`f * avg_wet_elev + (1 - f) * avg_land_elev`. Dividing by `f` -- the
pattern `isomip_plus`'s `TopoRemap` uses -- only recovers `avg_wet_elev`
if `avg_land_elev` is approximately zero, which holds for
`isomip_plus`'s idealized near-sea-level ice-shelf domain but not for
Great Lakes shorelines, where land elevation near the coast is tens of
meters positive. `GreatLakesTopoRemap._preprocess` therefore multiplies
`base_elevation` by the `(elevation < 0)` mask *before* the conservative
remap, zeroing land's contribution so the integral only accumulates over
the wet area; dividing the remapped result by the remapped `oceanFraction`
then exactly recovers the wet-area average (the same `var * mask`
pre-multiplication `polaris/tasks/e3sm/init/topo/remap/mask.py` uses for
its `*_masked_*` auxiliary fields, applied here to `base_elevation`
itself). Without this, a cell at 56% ocean fraction in Lake Erie --
barely above `min_ocean_fraction` -- showed a spurious `base_elevation`
of +26 m instead of a physically sensible shallow nearshore depth.

### Algorithm Design: Bathymetry coarsening before remap

Date last modified: 2026/10/23

`mbtempest`'s conservative-remap weight computation costs scale with
source point count, not target mesh size. NOAA's 3 arc-second (~90 m)
per-lake DEMs run into the tens of millions of points -- for Lake Erie,
remapping at native resolution onto the ~6,300-cell culled mesh ran over
55 minutes with no end in sight. `CoarsenBathymetryStep` block-averages
the DEM by `coarsen_factor` (default 3, i.e. ~90 m to ~270 m) via
`xarray.Dataset.coarsen(...).mean()` before the final
`GreatLakesTopoMap`/`GreatLakesTopoRemap` pair, cutting Erie's source
point count from 17.3M to 1.92M and the weight-computation step from
55+ minutes to under 5. The mesh step's cell-width shaping and
`WetMaskOnMeshStep`'s point sampling both use the native-resolution file
directly and are unaffected, since neither is the slow part. `mbtempest`
reports no completion fraction during weight computation (confirmed via
`mbtempest --verbose`, which only adds load/conversion timing and
per-event degenerate-triangle diagnostics, not a progress counter), so
`GreatLakesTopoMap` logs the source/target element counts once up front
and an elapsed-time heartbeat (`remap_progress_interval_s`) during the
run, as a liveness check rather than a true percentage.

## Implementation

### Implementation: Package layout

Date last modified: 2026/10/22

```
polaris/mesh/spherical/graded_grid.py          # build_graded_lon_lat()
polaris/tasks/ocean/great_lakes/
├── lakes.py            # GREAT_LAKES = [superior, michigan_huron, ...]
├── task.py             # GreatLakesMeshTask: wires one basin's step graph
├── great_lakes_mesh.cfg / great_lakes_topo.cfg
├── mesh/spherical.py   # GreatLakesSphericalMesh(QuasiUniformSphericalMeshStep)
├── topo/
│   ├── bathymetry.py   # PrepareBathymetryStep: decompress, standardize, merge
│   ├── mask.py         # WetMaskOnMeshStep: cheap pre-cull wet/dry mask
│   ├── coarsen.py      # CoarsenBathymetryStep: block-average before remap
│   ├── map.py          # GreatLakesTopoMap(MappingFileStep)
│   └── remap.py        # GreatLakesTopoRemap
└── cull.py             # GreatLakesCullMesh
```

`add_great_lakes_tasks(component)` (in `polaris/tasks/ocean/great_lakes/__init__.py`)
registers one `GreatLakesMeshTask` per basin, each with its own
`great_lakes_mesh.cfg`, at
`ocean/spherical/great_lakes/<basin>/mesh`.

### Implementation: Bathymetry input

Date last modified: 2026/10/23

`bathymetry_filenames` in `[great_lakes_topo]` is a plain, user-set,
comma- or whitespace-separated list of local paths (`<<<missing>>>` by
default), read at run time rather than setup time, following the pattern
`ComputeFeatureMasksStep` uses for its `mesh_filename` option -- this
avoids failing at `polaris setup` before the user has edited the work-dir
config. `PrepareBathymetryStep` transparently gunzips each source file if
needed (NOAA's Great Lakes DEMs ship as gzipped GMT NetCDF grids) and
renames its coordinate/elevation variables to `lon`/`lat`/`elevation`;
when more than one source is given, it merges them (see "Algorithm
Design: Merging Michigan and Huron's bathymetry" above). Most basins need
only one file; `michigan_huron` needs both `michigan_lld.grd.gz` and
`huron_lld.grd.gz`.

## Decisions

### Decision: Michigan and Huron combined into one basin

Lake Michigan and Lake Huron share the same water-surface elevation and
are connected by the Straits of Mackinac -- they are the same lake
hydrologically, unlike the other lake pairs. Unlike the remaining
connecting channels, their downloaded DEMs actually cover the connection
(see "Algorithm Design: Merging Michigan and Huron's bathymetry"), so
there was no reason to keep them as two artificially separate
sub-domains once that was confirmed. `michigan_huron` replaces the
`michigan`/`huron` entries in `GREAT_LAKES` as a single combined task.

### Decision: Four independent basins otherwise, no other connecting channels, no St. Lawrence

The remaining bathymetry DEMs used (NOAA's per-lake "lld" grids) do not
cover the other connecting channels (St. Marys River, St. Clair
River/Lake St. Clair/Detroit River, Niagara River), so Superior,
Michigan-Huron, Erie, and Ontario are disconnected rasters as provided.
Bridging them with idealized channel bathymetry, or sourcing real channel
surveys, and extending the domain to the St. Lawrence River, are
deferred. Polaris has no open-boundary or sponge-layer infrastructure
today, so a St. Lawrence connection will most likely be represented later
as a river-discharge flux at Lake Ontario's outlet (as E3SM's river model
couples to the ocean elsewhere) rather than a true open boundary.

### Decision: No per-lake water-level datum handling

NOAA's "lld" grids already have each lake's own IGLD85 water-level datum
subtracted (confirmed from the files' `history` attribute, e.g.
`grdmath erie_igld.grd 173.5 SUB = erie_lld.grd`), so `elevation < 0` is a
uniform wet/dry test across all five files. No per-lake threshold or
lake-ID raster is needed.

## Testing

### Testing and Validation: Graded grid

Date last modified: 2026/10/22

`tests/mesh/spherical/test_graded_grid.py` checks that
`build_graded_lon_lat` is monotonic and spans the globe, is uniform within
the fine region, grows geometrically outside it, and reduces to a uniform
grid when the fine and coarse resolutions match.

### Testing and Validation: End-to-end mesh

Date last modified: 2026/10/26

First run end-to-end for Lake Erie via `polaris setup` + `polaris serial`
on NERSC Perlmutter (`pm-cpu`, gnu/mpich), with the downloaded NOAA
`erie_lld.grd.gz` bathymetry, to validate the step graph and the
`coarsen_factor` performance fix (at the pre-fix `coastal_resolution_km` =
`slope_resolution_km` = 0.925 km): the full step graph (`bathymetry` ->
`base_mesh` -> `wet_mask` -> `cull_mesh` -> `coarsen` -> `topo_map_culled`
-> `topo_remap_culled`) completed successfully, producing a 6,289-cell
culled mesh with `base_elevation` ranging from -61.8 m to -0.09 m (all
negative, as expected after the topography normalization fix -- see
"Algorithm Design: Topography normalization") and `oceanFraction` > 0 for
6,287 of 6,289 cells (mean 0.97). With `coarsen_factor=3`, the
weight-computation step (`topo_map_culled`, 17.3M native DEM points
coarsened to 1.92M) completed in 4 m 47 s, versus 55+ minutes with no
completion in sight at native resolution before the coarsening step was
added. Ontario was also run at this point (full task runtime 2 m 16 s).

Superior was run end-to-end three times as the slope/trough fix was
developed, each still at `coastal_resolution_km` = `slope_resolution_km`
= 0.925 km:

1. With coastal-distance refinement only: deepest cell -350 m.
2. Adding slope-based refinement (`_fill_enclosed_steep_regions` not yet
   added): deepest cell -367.5 m at (46.91N, 86.60W), over a native-DEM
   pixel at -405.1 m, with an effective resolution there of only about
   2.83 km, since the trough floor's own slope never triggered refinement
   (see "Algorithm Design: Slope-based refinement").
3. Adding `_fill_enclosed_steep_regions` (`slope_fill_radius_km` = 2.0):
   mesh cell count rose from 19,553 to 20,634 (+5.5%), the trough floor's
   resolution reached 0.94 km (vs. 2.83 km before), and the deepest cell
   reached -389.7 m -- about 22 m deeper than before, and within 15.4 m
   of the native-DEM truth, though short of the idealized ~-399 m
   area-averaging estimate (the real JIGSAW cell isn't centered exactly on
   the single deepest native pixel the way that estimate assumed).

`coastal_resolution_km` and `slope_resolution_km` were then both dropped
to 0.8125 km (`max_resolution_km / 4`; see "Requirement: Mesh resolution").
Idealized area-averaging at Superior's trough gives about -401.4 m for
either 0.8125 km or 0.75 km (they round to the same averaging window on
the 3x-coarsened remap source, so there's no benefit to going finer than
0.8125 km without also reducing `coarsen_factor`).

All four basins were rerun end-to-end at `coastal_resolution_km` =
`slope_resolution_km` = 0.8125 km, `slope_fill_radius_km` = 2.0 (each on
an interactive or `shared`-QOS Perlmutter allocation, with
`MPICH_GPU_SUPPORT_ENABLED=0` -- the shell's default of 1 otherwise aborts
`mbconvert`, which isn't linked with GPU/GTL support, on a CPU-only node):

- **Ontario**: succeeded in 2 m 45 s; 5,632 cells; deepest cell -225.7 m
  at 3.10 km resolution (its deepest point is a broad basin floor, not an
  enclosed trough, so the slope/trough fix has no effect there).
- **Erie**: succeeded in 10 m 20 s; 6,936 cells (vs. 6,289 at 0.925 km,
  +10.3%, from the finer shoreline target, not the trough fix); deepest
  cell -61.8 m at 3.23 km resolution (also a broad floor, not a trough).
- **Superior**: succeeded in 42 m 39 s; 23,433 cells (vs. 20,634 at
  0.925 km with the trough fix, vs. 19,553 with coastal-distance
  refinement only -- a cumulative +19.8%); deepest cell -401.2 m at
  (46.9195N, 86.5955W), resolution 0.83 km -- within 4 m of the native-DEM
  truth (-405.1 m) and matching the -401.4 m idealized area-averaging
  estimate almost exactly, a much tighter match than the 0.925 km case's
  -389.7 m against its -399 m estimate.
- **Michigan-Huron**: succeeded in 1 h 23 m 40 s (longest of the four,
  from its larger combined domain -- `topo_map_culled` alone took
  1 h 10 m 24 s); 33,884 cells; deepest cell -271.4 m at (44.6897N,
  86.7345W), resolution 3.28 km (a broad basin floor near Lake Michigan's
  own deep basin, not an enclosed trough, so again outside what the
  slope/trough fix targets).

All four basins now build successfully end-to-end at the
`coastal_resolution_km` = `slope_resolution_km` = 0.8125 km,
`slope_fill_radius_km` = 2.0 configuration above. Each lake's single
globally-deepest point is a broad basin floor for Erie, Ontario, and
Michigan-Huron, and only Superior's is an enclosed steep-walled trough --
but checking only the deepest point per lake turned out to be incomplete:
Michigan-Huron has its own enclosed trough elsewhere (the Manitou
Passage, already well captured at these defaults) and several smaller,
gentler-walled basins in Georgian Bay that these shared defaults miss
entirely, and Ontario's own deepest point is under-refined for an
unrelated reason (`slope_fill_radius_km` falling just short of its
trough wall's gap). See "Algorithm Design: Topographic length-scale
refinement" for the diagnosis and fix: a second, scale-aware refinement
criterion run alongside the (unchanged) slope-based one above.

All four basins were rerun end-to-end again with that length-scale
criterion added (`length_scale_threshold_low_km = 1.5`,
`length_scale_threshold_high_km = 4.0`,
`length_scale_fill_radius_km = 3.0`; Superior and Michigan-Huron via a
`shared`-QOS Perlmutter batch job, account `m3780`; Erie and Ontario
interactively; all with `MPICH_GPU_SUPPORT_ENABLED=0`):

- **Michigan-Huron**: 33,884 -> 45,488 cells (+34.2% over the baseline
  above); the Georgian Bay pockets' bias falls from 24.3 m/5.4 m to
  2.8 m/2.2 m.
- **Ontario**: 5,632 -> 6,992 cells (+24.2%); deepest-point bias falls
  from 18.2 m to 5.7 m (-225.7 m to about -238.2 m, vs. -243.9 m truth).
- **Superior**: 23,433 -> 30,759 cells (+31.3%); deepest cell goes from
  -401.2 m to -401.4 m, essentially unchanged (truth -405.1 m) -- no
  regression, just added cost from sharing the same config as the lakes
  that need it.
- **Erie**: 6,936 -> 7,759 cells (+11.9%); its bathymetry has no real
  trough to fix, so this is pure added cost with no accuracy benefit --
  the accepted trade-off for keeping one config shared by all four
  basins rather than special-casing the three that need it.

See "Algorithm Design: Topographic length-scale refinement" for the full
derivation.
