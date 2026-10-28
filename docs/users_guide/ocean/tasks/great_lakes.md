(ocean-great-lakes)=

# Great Lakes

These tasks build a regional Omega/MPAS-Ocean mesh for a Laurentian Great
Lakes basin from a bathymetry DEM. They produce a mesh and bathymetry
only; there is no initial condition, forcing, or forward-run step yet.

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

## Example

```bash
polaris setup -t ocean/spherical/great_lakes/erie/mesh -w great_lakes_erie
```

Edit `great_lakes_erie/ocean/spherical/great_lakes/erie/mesh/great_lakes_mesh.cfg`
to point `bathymetry_filenames` at your downloaded DEM(s), then run the
task. Missing required options are reported when the task runs, after the
work-dir config has been edited.
