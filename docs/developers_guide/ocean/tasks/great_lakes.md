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
   multi-source basin.
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

Not yet implemented: bathymetry for the remaining connecting channels (St.
Marys River, St. Clair River/Lake St. Clair/Detroit River, Niagara River,
so the four basins cannot be joined into one mesh), initial conditions,
forcing, and a forward-run step.
