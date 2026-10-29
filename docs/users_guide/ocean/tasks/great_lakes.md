(ocean-great-lakes)=

# Great Lakes

These tasks build a regional Omega/MPAS-Ocean mesh for a Laurentian Great
Lakes basin from a bathymetry DEM, merge all four basins into a single
mesh, and build an idealized initial condition and a short forward run on
it.

Each basin is built and culled as its own closed sub-domain, except for
`michigan_huron`: Lake Michigan and Lake Huron are hydrologically a
single lake (same water-surface elevation, connected by the Straits of
Mackinac), and their downloaded DEMs tile together exactly, so they are
combined into one connected mesh. The remaining connecting channels (St.
Marys River, St. Clair River/Lake St. Clair/Detroit River, Niagara River)
and the St. Lawrence River are not yet included. See
`docs/design_docs/great_lakes_mesh.md` for why.

These tasks support both MPAS-Ocean and Omega.

(ocean-great-lakes-mesh)=

## mesh

### description

`ocean/spherical/great_lakes/<basin>/mesh` builds a culled MPAS mesh and
bathymetry for one basin, where `<basin>` is one of `superior`,
`michigan_huron`, `erie`, `ontario`.

### mesh

The mesh is spherical and quasi-uniform in the sense of
`polaris.mesh.spherical.QuasiUniformSphericalMeshStep`, with cell width
varying by distance to the lake's own bathymetry-derived shoreline rather
than by latitude: 0.8125 km at the shoreline, relaxing to a cap of 3.25 km
away from shore, in the `[great_lakes_mesh]` config section below.

### vertical grid

N/A (mesh-only task; no initial condition is produced yet).

### initial conditions

N/A

### forcing

N/A

### time step and run duration

N/A

### config options

```cfg
# config options for building a JIGSAW mesh for a single Great Lake
[great_lakes_mesh]

# cell width at the shoreline, km. 0.8125 = max_resolution_km / 4, and is
# close to 3x the ~270 m pixel spacing of the coarsened bathymetry used
# for the final remap
coastal_resolution_km = 0.8125

# cap / background cell width within the lake's own region, km
max_resolution_km = 3.25

# distance, km, over which cell width ramps from coastal_resolution_km to
# max_resolution_km
coastal_transition_km = 10.0

# exponent of the coastal transition blend; 1.0 is a linear ramp
coastal_transition_exponent = 1.0

# cell width, km, imposed far from the lake to keep JIGSAW tractable
background_width_km = 100.0

# buffer, degrees, added around the lake's own bathymetry extent
region_buffer_deg = 0.5

# distance, km, over which cell width relaxes to background_width_km
background_transition_km = 200.0

# resolution, km, applied wherever bathymetric slope exceeds
# slope_threshold_high (e.g. steep trench walls far from shore that the
# coastal-distance refinement above never reaches)
slope_resolution_km = 0.8125

# bathymetric slope (dimensionless: m of depth change per m of distance)
# below which no slope-based refinement is applied
slope_threshold_low = 0.02

# bathymetric slope above which slope_resolution_km fully applies
slope_threshold_high = 0.1

# exponent of the slope-based blend; 1.0 is a linear ramp
slope_transition_exponent = 1.0

# radius, km, of the morphological closing applied to the steep mask, so
# that low-slope trough floors fully enclosed by steep walls are refined
# to slope_resolution_km too, even though their own local slope never
# crosses slope_threshold_high
slope_fill_radius_km = 2.0

# resolution, km, applied wherever the topographic length scale (see
# length_scale_threshold_low_km below) drops to length_scale_threshold_low_km
# or less. Matches coastal_resolution_km/slope_resolution_km by default
length_scale_resolution_km = 0.8125

# bathymetric length scale L = he / |grad(H)| (he = harmonic-mean depth of
# two neighboring points, a scale-aware alternative to a raw slope value)
# at/below which length_scale_resolution_km fully applies
length_scale_threshold_low_km = 1.5

# length scale above which no length-scale-based refinement is applied
length_scale_threshold_high_km = 4.0

# exponent of the length-scale blend; 1.0 is a linear ramp
length_scale_transition_exponent = 1.0

# radius, km, of the morphological closing applied to the short-length-
# scale mask, analogous to slope_fill_radius_km above
length_scale_fill_radius_km = 3.0

# resolution, degrees, of the fine part of the graded JIGSAW sizing grid
fine_grid_resolution_deg = 0.0025

# resolution, degrees, of the coarse far-field part of the graded grid
coarse_grid_resolution_deg = 1.0

# per-step growth factor between the fine region and the coarse far field
grid_growth_factor = 1.15
```

The `slope_threshold_*`/`slope_fill_radius_km` defaults above were
calibrated against Lake Superior's and Lake Erie's native bathymetric
slope distributions, and on their own miss real troughs in lakes with
gentler overall relief (Lake Michigan-Huron's native DEM is roughly half
as steep as Superior's/Erie's at every percentile). Rather than
overriding the slope options per lake, the `length_scale_*` options above
add a second, scale-aware refinement criterion -- a topographic length
scale `L = he / |grad(H)|` -- that responds to *relative* depth change
instead of a raw slope value, so one shared set of defaults works for all
four basins without per-lake tuning. See "Algorithm Design: Topographic
length-scale refinement" in `docs/design_docs/great_lakes_mesh.md` for
the physical motivation and validated per-lake results.

```cfg
# config options for preparing and remapping a single Great Lake basin's
# bathymetry
[great_lakes_topo]

# path(s) to the basin's bathymetry file(s), comma- or whitespace-
# separated -- set this before running. Most basins need only one file;
# michigan_huron needs both michigan_lld.grd.gz and huron_lld.grd.gz
# (merged onto one grid; they tile together exactly at the Straits of
# Mackinac).
bathymetry_filenames = <<<missing>>>

# longitude, degrees, east of which native bathymetry is discarded before
# meshing (keeps lon <= lon_max_deg). Blank means no clipping. Set to
# -84.34 for `superior` only, to exclude the Whitefish Bay/St. Marys River
# approach to the Soo Locks at Sault Ste. Marie -- water hydraulically
# controlled by the locks, not part of the open lake; see "Decision: Lake
# Superior's domain stops at the Soo Locks" in
# `docs/design_docs/great_lakes_mesh.md`.
lon_max_deg =

# excludes a rectangular corner of the basin's bathymetry: source pixels
# with longitude < exclude_lon_max_deg AND latitude > exclude_lat_min_deg
# are set to land. Unlike lon_max_deg (which clips the whole grid at one
# longitude across every latitude), this only removes a corner, so it can
# isolate a connecting-channel approach that doesn't span the basin's
# full latitude range. Both must be set together; blank means no
# exclusion. Set to -83.96/46.28 for `michigan_huron` only, to exclude
# Huron's side of the same Soo Locks approach Superior's `lon_max_deg`
# excludes -- the two basins' independently built meshes were found to
# physically overlap there otherwise; see "Decision: Lake Huron's domain
# stops short of the Soo Locks too" in
# `docs/design_docs/great_lakes_mesh.md`.
exclude_lon_max_deg =
exclude_lat_min_deg =

# variable names in the source file(s)
elevation_var = z
lon_var = x
lat_var = y

# minimum ocean fraction of a remapped MPAS cell below which it is culled
min_ocean_fraction = 0.5

# factor by which the native-resolution bathymetry is block-averaged
# before the final conservative remap (mbtempest's weight-computation
# cost scales with source point count, so this keeps a tens-of-millions-
# point native DEM tractable while remaining several times finer than
# coastal_resolution_km). Set to 1 to disable and remap at native
# resolution.
coarsen_factor = 3

# how often, in seconds, the remap step logs a liveness heartbeat while
# mbtempest is computing weights (it reports no completion fraction, so
# this is only evidence the step hasn't hung, not a progress estimate)
remap_progress_interval_s = 60
```

### cores

The number of cores is the default for `polaris.remap.MappingFileStep`
(one task, `ntasks=1`); `coarsen_factor` keeps even the largest basin's
coarsened bathymetry grid small enough not to need parallel remapping.

(ocean-great-lakes-combined-mesh)=

## combined/mesh

### description

`ocean/spherical/great_lakes/combined/mesh` merges all four basins'
culled meshes into a single mesh spanning all four, and adds the Coriolis
parameter, vector-reconstruction weights, and a SCRIP file -- all
required before the mesh can be used as an Omega/MPAS-Ocean input, none
of which any single basin's mesh task produces. Depends on all four
`ocean/spherical/great_lakes/<basin>/mesh` tasks.

### mesh

A disjoint union of the four basins' culled meshes: cell/edge/vertex
arrays concatenated with renumbered IDs, no connectivity created between
basins (each remains its own closed sub-domain, as it already was). See
`docs/design_docs/great_lakes_mesh.md` for why this is correct rather
than a simplification.

### visualization

The `viz` step plots the combined mesh: full-domain cell width and
bathymetry (filled, with no cell outlines), plus four zoomed-in views
with crisp, individually visible cell outlines over bathymetry at the
Soo Locks, the Straits of Mackinac, Georgian Bay, and Lake Superior's
southeast glacial trough -- places where the combined mesh's
construction is easiest to check by eye, or (for the trough) where the
mesh's slope/length-scale refinement is best seen. Outlines are not
drawn on the full-domain plots: at 89,691 cells over the combined domain,
outlines at that scale blend into a solid wash rather than showing
individual cells. The zoomed views are colored by bathymetry rather than
cell width since depth is the more informative field at that scale; each
has its own colormap `vmin`, tuned to its own depth range (see the config
options below). `debug_summary.txt` records scalar diagnostics, including
each zoomed view's span-to-median-cell-width ratio as a legibility check.

```{image} images/great_lakes_combined_cell_width.png
:align: center
:width: 500px
```

```{image} images/great_lakes_combined_bathymetry.png
:align: center
:width: 500px
```

```{image} images/great_lakes_combined_mesh_outline_georgian_bay.png
:align: center
:width: 500px
```

```{image} images/great_lakes_combined_mesh_outline_superior_trough.png
:align: center
:width: 500px
```

### config options

```cfg
# config options for the mesh merging the four Great Lakes basins
[coriolis]

# Coriolis type for the merged mesh
type = spherical


# config options for visualizing the combined Great Lakes mesh
[great_lakes_combined_viz]

# the (lon_min, lon_max, lat_min, lat_max) extent, degrees, of the
# full-domain plots -- the union of the four basins' extents
extent = -92.4, -75.8, 41.2, 49.7


# colormap for the combined mesh's cell width.  Banded rather than
# continuous (12 steps of 0.25 km over the 0.5-3.5 km range below) so the
# coastal, slope, and interior resolution bands read as distinct bands
[viz_combined_cell_width]
colormap_name = viridis
colormap_levels = 12
norm_type = linear
norm_args = {'vmin': 0.5, 'vmax': 3.5}
colorbar_ticks = 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5
under_color = black
over_color = yellow


# colormap for the remapped bathymetry.  vmin is just below the combined
# mesh's deepest cell (about -401 m, in Lake Superior)
[viz_combined_bathymetry]
colormap_name = cmo.deep_r
norm_type = linear
norm_args = {'vmin': -410.0, 'vmax': 0.0}
under_color = black
over_color = tan

# colormaps for the four zoomed mesh-outline views, each with its own
# vmin tuned to that zoom's own depth range
[viz_combined_bathymetry_soo_locks]
colormap_name = cmo.deep_r
norm_type = linear
norm_args = {'vmin': -200.0, 'vmax': 0.0}
under_color = black
over_color = tan

[viz_combined_bathymetry_straits_of_mackinac]
colormap_name = cmo.deep_r
norm_type = linear
norm_args = {'vmin': -140.0, 'vmax': 0.0}
under_color = black
over_color = tan

[viz_combined_bathymetry_georgian_bay]
colormap_name = cmo.deep_r
norm_type = linear
norm_args = {'vmin': -200.0, 'vmax': 0.0}
under_color = black
over_color = tan

[viz_combined_bathymetry_superior_trough]
colormap_name = cmo.deep_r
norm_type = linear
norm_args = {'vmin': -410.0, 'vmax': 0.0}
under_color = black
over_color = tan
```

(ocean-great-lakes-combined-init)=

## combined/init

### description

`ocean/spherical/great_lakes/combined/init` builds a vertical grid and an
idealized initial condition on the combined mesh. Depends on
`ocean/spherical/great_lakes/combined/mesh`.

### vertical grid

A 40-level `tanh_dz` z-star grid, `bottom_depth = 420.0` m (with margin
over the merged mesh's true ~401 m deepest point), `partial_cell_type =
partial`. Cells shallower than `min_bottom_depth`, and cells with a
fill-value/NaN remapped elevation, are clamped to that floor before
`bottomDepth` is computed.

### initial conditions

Idealized and quiescent (zero velocity), not a climatology: uniform
near-surface temperature, a `tanh` thermocline, and a cold bottom layer,
with one background salinity value per basin. See
`docs/design_docs/great_lakes_mesh.md` for why, and for the known
shoreline artifact the depth clamp introduces.

### forcing

N/A -- no forcing is applied; see the `combined/forward` task below.

### config options

```cfg
# config options for the idealized Great Lakes initial condition, built on
# the combined mesh from polaris.tasks.ocean.great_lakes.combine

# Options related to the ocean component
[ocean]

# Equation of state type -- teos-10 is required, not a simplification;
# see docs/design_docs/great_lakes_mesh.md
eos_type = teos-10

[vertical_grid]
grid_type = tanh_dz
vert_levels = 40
bottom_depth = 420.0
min_layer_thickness = 1.0
max_layer_thickness = 20.0
coord_type = z-star
partial_cell_type = partial
min_pc_fraction = 0.1
min_vert_levels = 3

[great_lakes_init]

# positive-down depth, m, that bottomDepth is clamped to
min_bottom_depth = 5.0

# idealized summer-stratified rest state
surface_temperature = 20.0
bottom_temperature = 4.0
thermocline_depth = 15.0
thermocline_width = 5.0

# practical/absolute salinity, g/kg, one value per basin
salinity_superior = 0.05
salinity_michigan_huron = 0.1
salinity_erie = 0.2
salinity_ontario = 0.2
```

(ocean-great-lakes-combined-forward)=

## combined/forward

### description

`ocean/spherical/great_lakes/combined/forward` runs a short Omega/
MPAS-Ocean forward smoke test on the combined mesh and idealized initial
condition. Depends on `ocean/spherical/great_lakes/combined/init`. This
validates that the mesh/init pipeline actually runs; it is not a tuned or
physically realistic simulation.

### time step and run duration

`btr_dt`/`dt_btr_ratio`/`del2`/`del4` are derived from the combined
mesh's own finest edge and vertical grid, not Omega's global-mesh
defaults (which are unconditionally unstable at this resolution); see
`docs/design_docs/great_lakes_mesh.md`. The default run is a 6-hour smoke
test with hourly output.

### config options

```cfg
# config options for the Great Lakes combined-mesh forward smoke test
[great_lakes_forward]

# approximate number of cells in the combined mesh, used only to size the
# number of MPI tasks
approx_cell_count = 90998

# barotropic (external-gravity-wave) time step, s
btr_dt = 3.0

# ratio of the main (baroclinic) time step to btr_dt
dt_btr_ratio = 24

# Laplacian horizontal viscosity, m^2/s
del2 = 20.0

# Biharmonic horizontal viscosity, m^4/s
del4 = 1.0e6

# whether to disable the del4 (biharmonic) term entirely
disable_del4 = false

# implicit bottom drag coefficient (dimensionless)
bottom_drag_coeff = 1.0e-3

# length of the smoke-test run
run_duration_hours = 6

# how often output is written
output_interval_hours = 1
```

## Example

```bash
polaris setup -t ocean/spherical/great_lakes/erie/mesh -w great_lakes_erie
```

Edit `great_lakes_erie/ocean/spherical/great_lakes/erie/mesh/great_lakes_mesh.cfg`
to point `bathymetry_filenames` at your downloaded DEM(s), then run the
task. Missing required options are reported when the task runs, after the
work-dir config has been edited.

The combined mesh, init, and forward tasks need all four basins'
bathymetry files set up the same way:

```bash
polaris setup \
    -t ocean/spherical/great_lakes/superior/mesh \
       ocean/spherical/great_lakes/michigan_huron/mesh \
       ocean/spherical/great_lakes/erie/mesh \
       ocean/spherical/great_lakes/ontario/mesh \
       ocean/spherical/great_lakes/combined/mesh \
       ocean/spherical/great_lakes/combined/init \
       ocean/spherical/great_lakes/combined/forward \
    -w great_lakes_combined
```

Edit each basin's `great_lakes_mesh.cfg` as above, then run the tasks in
order (each depends on the previous).
