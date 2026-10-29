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

A second phase merges the four basins' culled meshes into one disjoint
mesh spanning all four, adds the Coriolis parameter and the
vector-reconstruction weights and SCRIP file a spherical mesh needs as an
ocean model input, and builds an idealized (quiescent, summer-stratified)
initial condition and a short Omega/MPAS-Ocean forward run on it -- closing
the gap this branch's own docs previously flagged ("mesh and bathymetry
only; no initial condition, forcing, or forward-run step yet").

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

### Requirement: Combined mesh

Date last modified: 2026/10/28

A single mesh spans all four basins, built by concatenating their culled
meshes with renumbered IDs and no connectivity created between them (each
basin stays a disjoint closed sub-domain, as it already is; see
"Decision: Four independent basins otherwise" above). The combined mesh
carries the Coriolis parameter, vector-reconstruction weights, and a
SCRIP file, none of which any single basin's culled mesh has, but which
`Ocean.write_horiz_mesh_dataset()` and the atm/ocean remapping that would
consume the SCRIP file both require.

### Requirement: Idealized initial condition and forward run

Date last modified: 2026/10/28

The combined mesh gets a vertical grid and an idealized, quiescent
(zero-velocity) initial condition -- a summer-stratified temperature
profile and a per-basin background salinity, not a real climatology, since
none was available (see "Decision: Idealized rest-state initial
condition" below) -- and a short Omega/MPAS-Ocean forward run that
exercises this state end to end. This closes the gap the mesh-only phase
above left open.

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

### Algorithm Design: Mesh merging

Date last modified: 2026/10/28

`mpas_tools.merge_grids.merge_grids()` already does almost exactly this
-- pairwise disjoint-mesh concatenation with per-variable ID-offset
renumbering -- and is the reference this step's offset table was taken
from. It cannot be called as-is, though: it pads the output's `maxEdges`
to the max of its two inputs, but then assigns each input's un-padded
connectivity array directly into the padded output, which raises a
shape-mismatch `ValueError` whenever the two inputs' `maxEdges` actually
differ. That mismatch is not hypothetical here: Superior and
Michigan-Huron's culled meshes use `maxEdges` = 9, Erie and Ontario's use
8, so every merge order hits this. `CombineGreatLakesMesh`
(`polaris/tasks/ocean/great_lakes/combine.py`) pads each basin's
`cellsOnCell`/`edgesOnCell`/`verticesOnCell`/`edgesOnEdge`/`weightsOnEdge`
to the combined mesh's own max `maxEdges` with the existing MPAS zero
("no neighbor") sentinel before concatenating, then applies the same
cell/edge/vertex ID offsets `merge_grids` does, generalized from two
inputs to four and implemented in xarray (`mpas_tools.io.open_dataset`/
`write_netcdf`) rather than raw `netCDF4.Dataset`, matching this
package's existing convention. `sort_mesh` is not rerun on the combined
mesh: each basin's own internal ordering is already sorted for
performance, and sorting a disjoint-component mesh would only reorder
which basin comes first, with no cross-basin locality to gain.

The merged mesh also gains the Coriolis parameter
(`polaris.coriolis.add_coriolis_to_dataset`), vector-reconstruction
weights (`polaris.mesh.reconstruct.compute_reconstruction_weights`, the
same function `polaris/tasks/e3sm/init/topo/cull/cull.py`'s
`CullMeshStep` already calls), and a SCRIP file. The SCRIP-writing part
of `CullMeshStep._create_scrip_file` (`MpasCellMeshDescriptor.to_scrip()`
plus an `ncks -5` conversion to NETCDF3_64BIT_DATA, since writing that
format directly is prohibitively slow) was factored out into
`polaris.mesh.scrip.write_mpas_scrip_file()` so both call sites share it
rather than duplicating the conversion step.

Checked directly against the current four basins' culled meshes
(2026/10/27 build; Superior's domain-clipping change in "Decision: Lake
Superior's domain stops at the Soo Locks" is included): 90,998 combined
cells (30,759 + 45,488 + 7,759 + 6,992), 289,891 edges, 198,583 vertices,
81 cells with a fill-value/NaN `base_elevation` (handled generically by
the initial-condition step below, not per basin), and a combined SCRIP
`grid_area` summing to about 243,000 km^2 -- matching the Great Lakes'
real total surface area, a useful end-to-end sanity check on the whole
merge-plus-SCRIP pipeline.

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
├── cull.py             # GreatLakesCullMesh
├── combine.py          # CombineGreatLakesMesh: merge + Coriolis/recon/SCRIP
├── init.py             # Init: idealized T/S/velocity initial condition
├── forward.py          # Forward: short Omega/MPAS-Ocean smoke-test run
├── great_lakes_combined.cfg / great_lakes_init.cfg / great_lakes_forward.cfg
└── forward.yaml
polaris/mesh/scrip.py   # write_mpas_scrip_file(), shared with e3sm/init/topo/cull
```

`add_great_lakes_tasks(component)` (in `polaris/tasks/ocean/great_lakes/__init__.py`)
registers one `GreatLakesMeshTask` per basin, each with its own
`great_lakes_mesh.cfg`, at `ocean/spherical/great_lakes/<basin>/mesh`, then
three more tasks that depend on all four: `GreatLakesCombinedMeshTask`
(`.../combined/mesh`), `GreatLakesCombinedInitTask` (`.../combined/init`),
and `GreatLakesCombinedForwardTask` (`.../combined/forward`). The per-basin
`GreatLakesMeshTask` exposes its `cull_mesh`/`topo_remap_culled` steps as
attributes so the combined task can depend on them without recreating
them (see "Implementation: Combined, init, and forward tasks" below).

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

### Implementation: Combined, init, and forward tasks

Date last modified: 2026/10/28

`CombineGreatLakesMesh` depends directly on each basin's `cull_mesh` and
`topo_remap_culled` steps (via the attributes noted above), not on a
shared/cached mesh database -- there is exactly one combined mesh, so
there is nothing to share. `Init` (`OceanIOStep`) follows
`polaris/tasks/ocean/realistic_global/init/initial_state.py`'s pattern
most closely of any existing task: it combines a mesh, external
bathymetry, and `reconstruction_weights.nc` the same way. `Forward`
(`OceanModelStep`) follows `polaris/tasks/ocean/single_column/forward.py`'s
pattern instead, reading `mesh.nc`/`vert_coord.nc`/`init.nc` as local
files from an upstream step rather than through the shared-mesh-database
mechanism `realistic_global`'s own forward step uses -- again, there is
only the one combined mesh.

`Forward.compute_cell_count()` returns a plain config number
(`[great_lakes_forward] approx_cell_count`) rather than a formula: every
other task's `compute_cell_count()` computes cell count from a resolution
parameter because those tasks cover a family of resolutions, but the
combined Great Lakes mesh is one fixed mesh, so a formula would be
circular. Time step, barotropic time step, and horizontal viscosity are
likewise plain `[great_lakes_forward]` config values, not scaled by a
`dt_per_km`-style formula (the convergence framework's own convention):
derived once, directly, against the combined mesh's own measured finest
edge (551.6 m) and the vertical grid's `bottom_depth` (420.0 m) --
`btr_dt` = 3.0 s against an ~8.6 s CFL limit, `del2` = 20.0 m^2/s against
an ~3800 s stability bound at the chosen 72 s main time step, `del4` =
1.0e6 m^4/s for roughly a 1-day damping time at the grid scale. Omega's
own global-mesh defaults (`ViscDel2` = 1.0e3, `ViscDel4` = 1.2e11) are
unconditionally unstable at this resolution and are not usable as-is;
`Forward` sets `update_eos=True` so `[ocean] eos_type` also governs the
run's actual `EosType`/`config_eos_type`, not only how the initial
condition's tracers were written.

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

### Decision: Lake Superior's domain stops at the Soo Locks

Date last modified: 2026/10/28

Unlike the other three lakes, Superior's downloaded NOAA DEM
(`superior_lld.grd.gz`) does not stop at its own shoreline: it extends
east to lon -84.0, into Whitefish Bay and the St. Marys River approach
toward Sault Ste. Marie -- water whose level and flow are hydraulically
controlled by the Soo Locks (and the adjacent International Bridge), not
representative of the open lake. This is the same class of feature the
previous decision already excludes for the other three connecting
channels (St. Marys River itself downstream of the locks, St. Clair
River/Lake St. Clair/Detroit River, Niagara River); Superior's DEM simply
happens to reach slightly past its own control point, where the other
lakes' DEMs don't reach theirs at all.

`PrepareBathymetryStep` (`polaris/tasks/ocean/great_lakes/topo/bathymetry.py`)
gained an optional `lon_max_deg` option in `[great_lakes_topo]`
(blank/unset by default, so the other three basins are unaffected):
when set, it clips the prepared `bathymetry.nc` to `lon <= lon_max_deg`
before anything downstream ever sees it, via a plain `ds.sel(lon=slice(None,
lon_max_deg))`. Since the mesh's bounding box
(`GreatLakesSphericalMesh.build_cell_width_lat_lon`), its shoreline
distance field, the pre-cull wet mask (`WetMaskOnMeshStep`), and the final
conservative remap all derive entirely from whatever `bathymetry.nc`
contains, clipping at this one upstream step is sufficient -- no other
code needed to change. `GreatLakesMeshTask` sets `lon_max_deg = -84.34`
(the user's estimate of the Soo Locks/Sault Ste. Marie bridge location)
only for `lake_name == 'superior'`, via a small `_LON_MAX_DEG_OVERRIDES`
dict in `task.py` -- a basin domain-extent decision, not a refinement
threshold, so it does not reopen the across-lake-consistency reasoning
behind the shared slope/length-scale refinement config (see "Algorithm
Design: Topographic length-scale refinement").

Checked directly against the real `bathymetry.nc`: this cutoff removes
103,594 of 14,305,485 native wet pixels (0.72%, about 606 km^2 of a
~83,700 km^2 lake). Checking wet-pixel longitude extent by latitude band
confirms the removed area is confined to lat 46.0-47.0, the Whitefish
Bay/St. Marys River corridor, not the open lake basin -- so this has no
effect on Superior's known trough or deepest point (see "Testing and
Validation: End-to-end mesh"), only a small reduction in domain size and
cell count.

### Decision: Lake Huron's domain stops short of the Soo Locks too

Date last modified: 2026/10/29

The decision above assumed "the other lakes' DEMs don't reach [their
control points] at all" -- true for Erie and Ontario, but not for
Michigan-Huron. While building `VizCombinedGreatLakesMesh`'s cell-width
plot, two basins that are supposed to be disjoint sub-domains (see
"Decision: Four independent basins otherwise") appeared visually
connected near Sault Ste. Marie. Direct analysis of the built
`combined_mesh.nc`, using each basin's actual cell polygons
(`verticesOnCell`/vertex coordinates) rather than just cell centers,
confirmed a real geometric overlap: 4 pairs of cells physically
overlapped (3 Superior cells, 2 Michigan-Huron cells), ~1.26 km^2 total,
confined to lon [-84.348, -84.329], lat [46.492, 46.508] -- right at
Superior's `lon_max_deg` cutoff. Huron's own downloaded DEM (the Huron
side, near the North Channel/St. Marys River mouth toward De Tour
Village) reaches into the same lock-controlled corridor Superior's does,
and its independently built mesh was never clipped to exclude it.

Unlike Superior's corridor, which spans its own full relevant latitude
range (so a plain `lon_max_deg` whole-grid clip is correct), Huron has
legitimate open-lake water -- the North Channel and Georgian Bay -- west
of Superior's cutoff longitude at lower latitudes. A whole-column clip
would remove that too, so `PrepareBathymetryStep` gained a second,
narrower mechanism: `exclude_lon_max_deg`/`exclude_lat_min_deg` in
`[great_lakes_topo]` (both blank/unset by default), which together mask
only a rectangular corner (`lon < exclude_lon_max_deg` AND
`lat > exclude_lat_min_deg`) to land via `.where()`, rather than dropping
whole rows/columns via `.sel()`. `GreatLakesMeshTask` sets
`exclude_lon_max_deg = -83.96` and `exclude_lat_min_deg = 46.28` only for
`lake_name == 'michigan_huron'`, via a `_MICHIGAN_HURON_NOTCH_EXCLUSION`
dict in `task.py`, mirroring `_LON_MAX_DEG_OVERRIDES`.

Checked directly against the real bathymetry: this exclusion removes
171,336 of 46,551,841 native wet pixels (0.37%). Depths in the excluded
corner are -8.7 to -0.3 m (median -3.8 m) -- shallow and river-like,
consistent with the lock-approach corridor rather than open lake, and
cells just south of the latitude cutoff have the same shallow profile,
confirming the cut falls within one homogeneous shallow feature rather
than across a lake/channel boundary. Rebuilding Michigan-Huron end to end
with the exclusion applied changes its culled-mesh cell count from 45,488
to 45,158 (-330 cells, -0.7%); the combined mesh (see "Algorithm Design:
Mesh merging") goes from 90,021 to 89,691 cells. `[great_lakes_forward]`'s
`approx_cell_count` was updated to match.

`CombineGreatLakesMesh` also gained `_check_no_basin_overlap`, run before
merging: for any two basins whose vertex-based bounding boxes intersect,
it builds real cell polygons and checks for true intersection with a
`shapely.STRtree`, raising a clear error naming the basins and
overlapping extent if any are found. The bounding-box pre-filter itself
had to use vertex coordinates rather than cell centers -- a first
implementation using `lonCell`/`latCell` bounds missed this overlap
entirely in a unit test, since a boundary cell's polygon can extend past
its own center by roughly a cell width, exactly the scale at which this
overlap occurred. Rerunning the check against the rebuilt meshes confirms
no overlap remains between any pair of basins.

### Decision: No per-lake water-level datum handling

NOAA's "lld" grids already have each lake's own IGLD85 water-level datum
subtracted (confirmed from the files' `history` attribute, e.g.
`grdmath erie_igld.grd 173.5 SUB = erie_lld.grd`), so `elevation < 0` is a
uniform wet/dry test across all five files. No per-lake threshold or
lake-ID raster is needed.

### Decision: Idealized rest-state initial condition, not a climatology

Date last modified: 2026/10/28

No Great Lakes temperature/salinity climatology was identified to build a
realistic initial condition from; building toward one is a much larger
follow-on effort, out of scope here. The initial condition is instead an
idealized, quiescent (zero-velocity) summer-stratified state: a uniform
near-surface temperature, a `tanh` thermocline, and a cold, well-mixed
bottom layer, all config options in `[great_lakes_init]`. This validates
the combined mesh, vertical grid, and model plumbing end to end; it is
not a tuned or physically realistic simulation, and should not be read as
one.

### Decision: TEOS-10, not a linear equation of state

Date last modified: 2026/10/28

`[ocean] eos_type = teos-10` is required, not a default left in place: a
linear equation of state has no pressure dependence, so it cannot
represent either the thermobaric effect (pressure's influence on the
thermal expansion coefficient, significant over Superior's ~400 m water
column) or the 4 degC freshwater density maximum that governs Great Lakes
stratification and seasonal overturn -- both central to why this task
exists on a GLS-capable Omega branch in the first place. This exercises
TEOS-10's IAPWS-95 pure-water formulation near its pure-water (near-zero
salinity) limit for the first time in this codebase; `eos_type = linear`
would be a safe debugging fallback if that limit misbehaves numerically,
but would also mean losing both effects, so it is not the default.
`Forward.__init__` passes `update_eos=True` so this option governs the
forward run's actual `EosType`/`config_eos_type`, not just how
`Init` writes the initial condition's tracers.

### Decision: Refine the vertical grid's near-surface resolution and layer growth rate

Date last modified: 2026/10/28

The original vertical grid (`vert_levels = 30`, `min_layer_thickness =
1.5`, `max_layer_thickness = 15.0`) grew too fast with depth: layer
thickness roughly tripled from layer 1 to layer 2 and again from layer 2
to layer 3, reaching its 15 m asymptote by layer 5. Only 5 of the 30
levels had their bottom shallower than 50 m, leaving the idealized
thermocline (`thermocline_depth = 15.0`, `thermocline_width = 5.0`)
sitting inside layers that were already 10-15 m thick, and the
`min_bottom_depth` floor (below) carried only 2 valid levels instead of
the intended 3 (see the superseded paragraph this replaces).

`vert_levels`/`min_layer_thickness`/`max_layer_thickness` were changed to
40/1.0/20.0 -- kept in the 30-40 range (rather than, e.g., 50, which
would resolve the thermocline band more finely still) to limit the
increase in column count to ~33% over the previous 30. This caps the
growth in thickness between adjacent layers at ~13% (down from up to
~2.7x), keeps the surface layer at 1.0 m, and keeps layer thickness
under ~2 m down to ~10 m depth, ~2.7 m at the 15 m thermocline center,
and ~7 m by 50 m depth (versus ~15 m by 50 m with the previous grid; 16
of 40 levels now have their bottom shallower than 50 m).

### Decision: `min_bottom_depth` = 5.0 m, sized to the vertical grid's own floor

Date last modified: 2026/10/28

Cells shallower than `[great_lakes_init] min_bottom_depth`, and the 81
cells (of 90,998) with a fill-value/NaN `base_elevation`, are clamped to
that floor before `bottomDepth` is computed. The chosen value, 5.0 m, is
sized to the vertical grid's own `min_vert_levels * min_layer_thickness`
(now 3.0 m, after the refinement above; 4.5 m previously), not deeper: a
deeper floor only enlarges the shoreline artifact this clamp introduces
(it affects 9.0% of cells at 5.0 m, versus 21.5% at a 10 m floor --
these two figures are a function of the remapped bathymetry only and are
unaffected by the vertical-grid refinement above) without buying
`init_vertical_coord` anything.

This does not fully eliminate the artifact `partial_cell_type = partial`
introduces on top of it. Under the previous grid, 12,181 of 90,998 cells
(13.4%) ended up with 2 valid vertical levels instead of 3 once the
bottom-most partial layer was trimmed under `min_pc_fraction`, though
none ended up with fewer than 2 or with zero -- i.e. no cell lost
vertical structure entirely, just somewhat more than the nominal minimum
intended. Replaying `compute_min_max_level_cell` and
`alter_bottom_depth` against the refined 1D grid at `bottomDepth =
min_bottom_depth = 5.0` m gives 5 valid levels rather than 3, a wider
margin than before, which should substantially shrink or eliminate that
shortfall -- **this has not yet been remeasured against the real merged
mesh's bathymetry distribution; do so, and update this paragraph's
percentage/cell-count, the next time the init step is run end-to-end.**

### Decision: Per-basin background salinity

Date last modified: 2026/10/28

The idealized salinity field is uniform with depth and within each
basin, but differs between basins (`[great_lakes_init]
salinity_<basin>`): 0.05 g/kg for Superior, 0.1 g/kg for Michigan-Huron,
0.2 g/kg for Erie and Ontario, reflecting the real Great Lakes' actual
range of background mineralization rather than treating the whole
combined mesh as one uniform value. `CombineGreatLakesMesh` adds a
`basinIndex` field (and a `basin_names` global attribute giving the
order) to `combined_topography.nc` precisely so `Init` can look up each
cell's basin; this is the only reason that field exists; it is not
needed by, and is not read by, the mesh-merging logic itself.

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

### Testing and Validation: Combined mesh, init, and forward

Date last modified: 2026/10/28

`tests/ocean/great_lakes/test_combine.py` unit-tests the merge logic
(`_merge_meshes`/`_pad_ragged_dims`/`_apply_id_offsets`) against small,
hand-built two- and three-basin fixtures, including one with mismatched
`maxEdges` (4 vs. 3) specifically to exercise the padding path that
`mpas_tools.merge_grids.merge_grids()` gets wrong (see "Algorithm Design:
Mesh merging"). It checks: combined cell/edge/vertex counts equal the
sum of inputs; no cell references an index outside its own basin's
offset range (no accidental cross-basin connectivity); padded columns
are the zero sentinel, not garbage; plain (non-connectivity) fields are
concatenated unchanged; and global attributes that agree between inputs
are carried forward.

`CombineGreatLakesMesh`, `Init`, and `Forward` were also each run
end-to-end via `polaris setup` + `polaris serial` on Perlmutter
(`pm-cpu`), seeded with the four basins' real culled meshes and remapped
bathymetry (the 2026/10/27 build referenced above):

- `combine`: produced a 90,998-cell `combined_mesh.nc` with `fCell`/
  `fEdge`/`fVertex` present, a `reconstruction_weights.nc` with no NaNs,
  and a `combined_mesh.scrip.nc` whose `grid_area` sums to the ~243,000
  km^2 sanity check noted above.
- `init`: produced `mesh.nc`/`vert_coord.nc`/`init.nc` with every
  variable `Ocean.write_horiz_mesh_dataset`/`write_vert_coord_dataset`/
  `write_initial_state_dataset` require present (`FCell`,
  `ReconWeightsCell`; `MinLayerCell`, `MaxLayerCell`, `BottomGeomDepth`,
  `VertCoordMovementWeights`, `RefPseudoThickness`; `Temperature`,
  `Salinity`, `NormalVelocity`, `PseudoThickness`, `SurfacePressure`).
  `Temperature` ranges 4.0-19.95 degC as configured; `Salinity` is
  exactly uniform within each basin at the configured per-basin value
  (0.05/0.1/0.2/0.2 g/kg); no cell ends up with zero valid vertical
  levels (see "Decision: `min_bottom_depth` = 5.0 m, sized to the
  vertical grid's own floor" for the 2-level shortfall this does not
  fully avoid). **This run predates the "Refine the vertical grid's
  near-surface resolution and layer growth rate" decision above and
  should be redone against the current 40-level grid** to confirm
  `Temperature`'s range is unaffected (it should be, since the
  vertical-grid change only affects how finely the same idealized
  profile is sampled) and to remeasure the 2-level shortfall.
- `forward`: `polaris setup` resolved `omega.yml` with the intended
  numeric values (`TimeStep` 0000_00:01:12, `BtrTimeStep`
  0000_00:00:03, `ViscDel2` 20, `ViscDel4` 1e+06, `BottomDragCoeff`
  0.001, `EosType: teos-10`, convective/shear mixing enabled, `StopCriterion`
  0000_06:00:00) and a dynamically-sized task count (456 target/44
  minimum, from `compute_cell_count()`'s 90,998). Actually executing the
  compiled model was not done in this session (no Omega build was
  available); this validates the full configuration pipeline up to that
  point, not the run's numerical behavior.
