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

Success is a culled MPAS mesh per basin with cell width 925 m at the
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
grid fine enough for 925 m resolution would have on the order of 10^8
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

Date last modified: 2026/10/23

Run end-to-end for Lake Erie via `polaris setup` +
`polaris serial` on NERSC Perlmutter (`pm-cpu`, gnu/mpich), with the
downloaded NOAA `erie_lld.grd.gz` bathymetry. The full step graph
(`bathymetry` -> `base_mesh` -> `wet_mask` -> `cull_mesh` -> `coarsen` ->
`topo_map_culled` -> `topo_remap_culled`) completed successfully,
producing a 6,289-cell culled mesh with `base_elevation` ranging from
-61.8 m to -0.09 m (all negative, as expected after the topography
normalization fix -- see "Algorithm Design: Topography normalization")
and `oceanFraction` > 0 for 6,287 of 6,289 cells (mean 0.97). With
`coarsen_factor=3`, the weight-computation step (`topo_map_culled`,
17.3M native DEM points coarsened to 1.92M) completed in 4 m 47 s, versus
55+ minutes with no completion in sight at native resolution before the
coarsening step was added.

Also run end-to-end for Ontario (full task runtime 2 m 16 s). The
`michigan_huron` combined basin is a substantially larger domain (lon
-88.0 to -79.68, lat 41.62 to 46.5, versus Erie's -83.7 to -78.5,
41.2-43.0) and is running as of this writing.

Superior has been run end-to-end three times as this fix was developed:

1. With coastal-distance refinement only: deepest cell -350 m.
2. Adding slope-based refinement (`_fill_enclosed_steep_regions` not yet
   added): deepest cell -367.5 m at (46.91N, 86.60W), over a native-DEM
   pixel at -405.1 m, with an effective resolution there of only about
   2.83 km, since the trough floor's own slope never triggered refinement
   (see "Algorithm Design: Slope-based refinement").
3. Adding `_fill_enclosed_steep_regions` (at `slope_resolution_km` =
   `coastal_resolution_km` = 0.925 km, `slope_fill_radius_km` = 2.0):
   mesh cell count rose from 19,553 to 20,634 (+5.5%), the trough floor's
   resolution reached 0.94 km (vs. 2.83 km before), and the deepest cell
   reached -389.7 m -- about 22 m deeper than before, and within 15.4 m
   of the native-DEM truth, though short of the idealized ~-399 m
   area-averaging estimate (the real JIGSAW cell isn't centered exactly on
   the single deepest native pixel the way that estimate assumed).

`coastal_resolution_km` and `slope_resolution_km` were then both dropped
to 0.8125 km (`max_resolution_km / 4`; see "Requirement: Mesh resolution").
Idealized area-averaging at this trough gives about -401.4 m for either
0.8125 km or 0.75 km (they round to the same averaging window on the
3x-coarsened remap source, so there's no benefit to going finer than
0.8125 km without also reducing `coarsen_factor`). Rerunning Superior
end-to-end with 0.8125 km to confirm the real-world result is still
pending.
