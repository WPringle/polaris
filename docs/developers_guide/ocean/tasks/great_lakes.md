(dev-ocean-great-lakes)=

# Great Lakes

{py:mod}`polaris.tasks.ocean.great_lakes` builds a culled, bathymetry-
informed mesh for a Laurentian Great Lakes basin. The design rationale,
including why `michigan_huron` combines two lakes while the other basins
are separate closed sub-domains, is in
`docs/design_docs/great_lakes_mesh.md`.

{py:func}`polaris.tasks.ocean.great_lakes.add_great_lakes_tasks` registers
one {py:class}`polaris.tasks.ocean.great_lakes.task.GreatLakesMeshTask` per
name in {py:data}`polaris.tasks.ocean.great_lakes.lakes.GREAT_LAKES`, at
`ocean/spherical/great_lakes/<basin>/mesh`. Each task wires the same step
graph:

1. {py:class}`polaris.tasks.ocean.great_lakes.topo.bathymetry.PrepareBathymetryStep`
   decompresses each user-supplied DEM if gzipped and renames its
   coordinate/elevation variables to `lon`/`lat`/`elevation`. The source
   path(s) (`bathymetry_filenames` in `[great_lakes_topo]`, comma- or
   whitespace-separated) are read at run time, not setup time, the same
   way
   {py:class}`polaris.tasks.mesh.spherical.feature_masks.ComputeFeatureMasksStep`
   handles its `mesh_filename` option, since it defaults to
   `<<<missing>>>` until the user edits the work-dir config. When more
   than one source is given (`michigan_huron`), they are merged onto one
   grid sized to the union of their extents -- see the design doc's
   "Merging Michigan and Huron's bathymetry" section. Every step after
   this one is basin-agnostic and needed no changes to support a
   multi-source basin. Two clipping mechanisms remove bathymetry that
   isn't part of the open lake: `lon_max_deg` drops whole longitude
   columns across every latitude (used for Superior's Soo Locks
   approach, which spans its own full relevant latitude range), while
   `exclude_lon_max_deg`/`exclude_lat_min_deg` mask only a rectangular
   corner via `.where()` rather than `.sel()` (used for Huron's side of
   the same corridor, since Huron has legitimate open-lake water --
   North Channel/Georgian Bay -- west of that longitude at lower
   latitudes, which a whole-column clip would also remove). The
   corner exclusion was added after the combined mesh's two
   independently built basins were found to physically overlap right at
   the Soo Locks; see the design doc's "Decision: Lake Huron's domain
   stops short of the Soo Locks too" section.
2. {py:class}`polaris.tasks.ocean.great_lakes.mesh.spherical.GreatLakesSphericalMesh`,
   a {py:class}`polaris.mesh.spherical.QuasiUniformSphericalMeshStep`
   subclass, builds the cell-width field. See
   {py:mod}`polaris.mesh.spherical.graded_grid` for the non-uniform JIGSAW
   sizing grid this depends on, and the design doc for the coastal-ramp
   and background-coarsening algorithm. A second, independent refinement
   based on bathymetric slope (`_bathymetric_slope`/`_slope_cell_width`)
   is combined with the coastal field via `np.minimum`, so steep features
   far from shore (e.g. Lake Superior's deep trenches, which a
   coastal-distance-only field left smoothed from -405 m to -350 m) get
   refined too. Since per-cell slope refinement alone still misses a flat
   trough floor enclosed by steep walls (slope-refined only up to -367.5 m
   for Superior), `_fill_enclosed_steep_regions` morphologically closes
   the steep-cell mask so such enclosed floors are refined to
   `slope_resolution_km` too, without expanding refinement along ordinary
   unenclosed steep terrain -- see the design doc's "Slope-based
   refinement" section.
3. {py:class}`polaris.tasks.ocean.great_lakes.topo.mask.WetMaskOnMeshStep`
   builds a cheap, point-sampled wet/dry mask on the uncull base mesh
   (nearest-neighbor lookup against the native DEM, gated to the DEM's
   extent, plus connectivity cleanup via
   {py:mod}`polaris.mesh.connectivity`) -- not a full conservative remap,
   since the cull decision doesn't need an area-accurate ocean fraction
   and the base mesh is mostly discarded by culling anyway.
4. {py:class}`polaris.tasks.ocean.great_lakes.cull.GreatLakesCullMesh`
   culls cells below `min_ocean_fraction`, reusing the
   `_land_mask_from_topo` idiom from
   {py:mod}`polaris.tasks.ocean.isomip_plus.mesh.cull` directly.
5. {py:class}`polaris.tasks.ocean.great_lakes.topo.coarsen.CoarsenBathymetryStep`
   block-averages the native DEM by `coarsen_factor` (`[great_lakes_topo]`,
   default 3) before the final conservative remap, since `mbtempest`'s
   weight-computation cost scales with source point count and the native
   DEM can run into the tens of millions of points.
6. {py:class}`polaris.tasks.ocean.great_lakes.topo.map.GreatLakesTopoMap`
   and {py:class}`polaris.tasks.ocean.great_lakes.topo.remap.GreatLakesTopoRemap`
   run the full conservative remap of the coarsened bathymetry onto the
   culled mesh, producing the final `base_elevation`/`oceanFraction`.
   `GreatLakesTopoMap` logs the source/target element counts up front and
   an elapsed-time heartbeat during weight computation
   (`remap_progress_interval_s`), since `mbtempest` reports no completion
   fraction. `GreatLakesTopoRemap` zeroes land's contribution to
   `base_elevation` (multiplying by the `(elevation < 0)` mask) before the
   remap, so a partially-land cell's final `base_elevation` is the
   area-weighted average depth of its wet portion only -- unlike
   `isomip_plus`'s `TopoRemap`, which divides the raw (land-included)
   remapped elevation by `oceanFraction` directly, an approximation only
   valid when land elevation is near zero. See the design doc's
   "Topography normalization" section for the derivation.

All steps in a task share one {py:class}`polaris.config.PolarisConfigParser`,
built from `great_lakes_mesh.cfg` and `great_lakes_topo.cfg` in
{py:mod}`polaris.tasks.ocean.great_lakes`, plus `polaris.mesh.spherical`'s
`spherical.cfg` and `polaris.remap`'s `mapping.cfg`.

## Combined mesh, init, and forward tasks

{py:func}`polaris.tasks.ocean.great_lakes.add_great_lakes_tasks` also
registers three tasks that depend on all four basins'
{py:class}`~polaris.tasks.ocean.great_lakes.task.GreatLakesMeshTask`
instances, which expose their `cull_mesh` and `topo_remap_culled` steps
as attributes for exactly this purpose:

1. {py:class}`polaris.tasks.ocean.great_lakes.task.GreatLakesCombinedMeshTask`
   (`ocean/spherical/great_lakes/combined/mesh`) runs
   {py:class}`polaris.tasks.ocean.great_lakes.combine.CombineGreatLakesMesh`,
   which merges the four basins' `culled_mesh.nc`/`topography_remapped.nc`
   into `combined_mesh.nc`/`combined_topography.nc` (`_merge_meshes`,
   `_pad_ragged_dims`, `_apply_id_offsets`; see the design doc's "Mesh
   merging" section for why `mpas_tools.merge_grids.merge_grids()` isn't
   called directly). Before merging, `_check_no_basin_overlap` verifies
   the disjoint-sub-domain assumption each basin's docstring asserts,
   rather than just trusting it: for any two basins whose bounding boxes
   intersect (built from vertex, not cell-center, coordinates, since a
   boundary cell's polygon can extend past its own center by roughly a
   cell width -- exactly the scale of the Soo Locks overlap this check
   was added to catch), it builds real cell polygons from
   `verticesOnCell` and checks for true intersection with a
   `shapely.STRtree`, raising a `ValueError` naming the basins and
   overlapping extent if any are found. `CombineGreatLakesMesh.run()`
   then regenerates `combined_graph.info` via
   {py:func}`polaris.model_step.make_graph_file`, adds the Coriolis
   parameter via {py:func}`polaris.coriolis.add_coriolis_to_dataset`, and
   writes `reconstruction_weights.nc` and `combined_mesh.scrip.nc` via
   {py:func}`polaris.mesh.reconstruct.compute_reconstruction_weights` and
   the new {py:func}`polaris.mesh.scrip.write_mpas_scrip_file` (factored
   out of {py:class}`polaris.tasks.e3sm.init.topo.cull.cull.CullMeshStep`,
   which now calls it too). It also adds a `basinIndex` field and a
   `basin_names` global attribute to `combined_topography.nc`, so the init
   step below can assign a per-basin background salinity. The task's
   second step, {py:class}`polaris.tasks.ocean.great_lakes.viz.VizCombinedGreatLakesMesh`,
   plots the merged mesh's cell width (mean of `dcEdge` over each cell's
   edges via `_cell_width_km`, not the `areaCell`-based hexagon estimate
   used elsewhere, since that estimate assumes a regular hexagon and
   understates width on the clipped coastal cells this plot exists to
   show) and bathymetry, with no cell outlines drawn on either. A
   full-domain outline plot isn't produced: at 90,021 cells, ink from the
   outlines would dominate the small coastal cells and wash into a solid
   color at any dpi, since the stroke width and cell size both scale with
   it. Instead, `MESH_OUTLINE_ZOOMS` draws three small boxes (the Soo
   Locks truncation, the Michigan/Huron bathymetry seam at the Straits of
   Mackinac, and the topographic-length-scale-refined troughs in Georgian
   Bay) with fully legible per-cell outlines, each one chosen to show a
   feature specific to how the combined mesh was built (see
   `docs/design_docs/great_lakes_mesh.md`). `debug_summary.txt` records
   each zoom's span-to-median-cell-width ratio so future boxes can be
   checked against the same rule of thumb that ruled out a full-domain
   version.
2. {py:class}`polaris.tasks.ocean.great_lakes.task.GreatLakesCombinedInitTask`
   (`ocean/spherical/great_lakes/combined/init`) runs
   {py:class}`polaris.tasks.ocean.great_lakes.init.Init`
   ({py:class}`~polaris.ocean.model.OceanIOStep`, following
   {py:class}`polaris.tasks.ocean.realistic_global.init.initial_state.InitialStateStep`'s
   pattern most closely of any existing step): clamps `bottomDepth` to
   `min_bottom_depth` (`_bottom_depth_from_elevation`), calls
   {py:func}`polaris.ocean.vertical.init_vertical_coord`, builds the
   idealized temperature profile (`_idealized_temperature`) and per-basin
   salinity (`_salinity_by_basin`, keyed by `basinIndex`), and adds zero
   velocity via
   {py:func}`polaris.ocean.init_state.add_quiescent_normal_velocity`
   before writing `mesh.nc`/`vert_coord.nc`/`init.nc`. T/S must be set on
   the dataset before `write_vert_coord_dataset`:
   {py:func}`polaris.ocean.vertical.diagnostics.pseudothickness_from_ds`
   silently returns `None` without them, which would leave
   `RefPseudoThickness` missing.
3. {py:class}`polaris.tasks.ocean.great_lakes.task.GreatLakesCombinedForwardTask`
   (`ocean/spherical/great_lakes/combined/forward`) runs
   {py:class}`polaris.tasks.ocean.great_lakes.forward.Forward`
   ({py:class}`~polaris.ocean.model.OceanModelStep`, following
   {py:class}`polaris.tasks.ocean.single_column.forward.Forward`'s
   pattern of reading the upstream step's local
   `mesh.nc`/`vert_coord.nc`/`init.nc` directly, rather than
   `realistic_global`'s shared-mesh-database mechanism, since there is
   only the one combined mesh). `compute_cell_count()` returns a plain
   `[great_lakes_forward] approx_cell_count` rather than a resolution
   formula, since the combined mesh is fixed, not a family of
   resolutions. `dynamic_model_config` templates `forward.yaml` with
   `btr_dt`/`dt_btr_ratio`/`del2`/`del4` derived directly against the
   combined mesh's own finest edge and vertical grid (see the design
   doc), and `update_eos=True` is passed to `super().__init__` so
   `[ocean] eos_type` also governs the run's `EosType`/`config_eos_type`.
   `forward.yaml`'s `Omega.IOStreams.InitialVertCoord.Filename` must be
   `vert_coord.nc`, matching `single_column`'s pattern (not
   `realistic_global`'s `init.nc`, which only works there because that
   task's `init.nc` holds the vertical-coordinate fields too): this task
   writes them to a separate `vert_coord.nc`
   ({py:class}`polaris.tasks.ocean.great_lakes.init.Init`), and pointing
   `InitialVertCoord` at `init.nc` instead crashed Omega with `Error
   reading BottomGeomDepth from InitialVertCoord` before any time steps
   ran.

Each of these three tasks builds its own
{py:class}`~polaris.config.PolarisConfigParser` (`great_lakes_combined.cfg`,
`great_lakes_init.cfg`, `great_lakes_forward.cfg` respectively, all in
{py:mod}`polaris.tasks.ocean.great_lakes`) rather than sharing the
per-basin mesh tasks' config, since none of their options are per-basin.
`GreatLakesCombinedForwardTask` also loads `great_lakes_init.cfg` for its
`[ocean] eos_type`, which must match the one the initial condition's
tracers were written with.

Not yet implemented: bathymetry for the remaining connecting channels (St.
Marys River, St. Clair River/Lake St. Clair/Detroit River, Niagara River,
so the four basins cannot be hydraulically connected), a realistic
temperature/salinity climatology (the initial condition above is
idealized), and surface forcing.
