#: The four basins built as Great Lakes meshes.
#:
#: Lake Michigan and Lake Huron are hydrologically a single lake: they
#: share the same water-surface elevation and are connected by the ~5 km
#: wide Straits of Mackinac. NOAA ships their bathymetry as two separate
#: per-lake grids, but as downloaded the two tile together exactly (same
#: 3 arc-second resolution, an identical shared boundary column right at
#: the strait), so they are combined into one ``michigan_huron`` task by
#: ``PrepareBathymetryStep`` (see
#: :py:mod:`polaris.tasks.ocean.great_lakes.topo.bathymetry`) rather than
#: built as two separate closed sub-domains.
#:
#: The remaining connecting channels (St. Marys River, St. Clair
#: River/Lake St. Clair/Detroit River, Niagara River) are not covered by
#: any of the downloaded bathymetry, so Superior, Michigan-Huron, Erie,
#: and Ontario remain their own closed sub-domains for now. Joining them
#: into one connected Great Lakes mesh is future work once
#: connecting-channel bathymetry is available.
GREAT_LAKES = ['superior', 'michigan_huron', 'erie', 'ontario']
