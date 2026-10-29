import numpy as np
import xarray as xr
from mpas_tools.io import open_dataset

from polaris import Step
from polaris.viz import plot_global_mpas_field

#: Zoom boxes (lon_min, lon_max, lat_min, lat_max), degrees, drawn with
#: visible cell outlines.  Each one is small enough that the outlines
#: stay legible (see ``_ZOOM_PAD_DEG`` and the ratio written to
#: ``debug_summary.txt``) and each shows a feature specific to the
#: combined mesh: the Soo Locks truncation of Superior's domain, the
#: seam where Michigan and Huron's separate bathymetry sources join, and
#: the topographic-length-scale refinement in Georgian Bay (see
#: docs/design_docs/great_lakes_mesh.md).
MESH_OUTLINE_ZOOMS = {
    'soo_locks': (-85.2, -83.8, 45.8, 47.0),
    'straits_of_mackinac': (-85.6, -84.2, 45.4, 46.3),
    'georgian_bay': (-82.0, -79.6, 44.4, 46.3),
}

#: Padding, in degrees, applied when selecting cells for a zoom box so
#: that cells right at the plotted edge are not clipped
_ZOOM_PAD_DEG = 0.05


class VizCombinedGreatLakesMesh(Step):
    """
    A step for visualizing the full combined Great Lakes mesh: cell
    width, bathymetry, and (in a few small regions where the outlines
    stay legible) the individual mesh cells.

    Attributes
    ----------
    combine_step : polaris.tasks.ocean.great_lakes.combine.CombineGreatLakesMesh
        The step that builds the combined mesh this step visualizes
    """  # noqa: E501

    def __init__(self, component, subdir, config, combine_step):
        """
        Create the step

        Parameters
        ----------
        component : polaris.Component
            The component the step belongs to

        subdir : str
            The subdirectory for the step

        config : polaris.config.PolarisConfigParser
            A shared config parser

        combine_step : polaris.tasks.ocean.great_lakes.combine.CombineGreatLakesMesh
            The step that builds the combined mesh
        """  # noqa: E501
        super().__init__(component=component, name='viz', subdir=subdir)
        self.set_shared_config(config, link='great_lakes_combined.cfg')
        self.combine_step = combine_step

        self.add_input_file(
            filename='combined_mesh.nc',
            work_dir_target=f'{combine_step.path}/combined_mesh.nc',
        )
        self.add_input_file(
            filename='combined_topography.nc',
            work_dir_target=f'{combine_step.path}/combined_topography.nc',
        )

        self.output_filenames = [
            'combined_cell_width.png',
            'combined_bathymetry.png',
            *(
                f'combined_mesh_outline_{name}.png'
                for name in MESH_OUTLINE_ZOOMS
            ),
            'debug_summary.txt',
        ]
        for filename in self.output_filenames:
            self.add_output_file(filename=filename)

    def run(self):
        """
        Run this step
        """
        config = self.config
        section = config['great_lakes_combined_viz']
        extent = tuple(section.getnumpy('extent'))

        with (
            open_dataset('combined_mesh.nc') as ds_mesh,
            open_dataset('combined_topography.nc') as ds_topo,
        ):
            width_km = _cell_width_km(ds_mesh)
            da_width = xr.DataArray(width_km, dims='nCells')
            base_elevation = ds_topo['base_elevation']

            descriptor = plot_global_mpas_field(
                mesh_filename='combined_mesh.nc',
                da=da_width,
                out_filename='combined_cell_width.png',
                config=config,
                colormap_section='viz_combined_cell_width',
                title='Combined Great Lakes mesh: cell width (km)',
                colorbar_label='cell width (km)',
                plot_land=False,
                extent=extent,
            )
            plot_global_mpas_field(
                mesh_filename='combined_mesh.nc',
                da=base_elevation,
                out_filename='combined_bathymetry.png',
                config=config,
                colormap_section='viz_combined_bathymetry',
                title='Combined Great Lakes mesh: base_elevation (m)',
                colorbar_label='elevation (m)',
                plot_land=False,
                extent=extent,
                descriptor=descriptor,
            )

            zoom_lines = []
            for name, box in MESH_OUTLINE_ZOOMS.items():
                cell_indices = _cells_in_box(
                    ds_mesh, _pad_box(box, _ZOOM_PAD_DEG)
                )
                plot_global_mpas_field(
                    mesh_filename='combined_mesh.nc',
                    da=da_width[cell_indices],
                    out_filename=f'combined_mesh_outline_{name}.png',
                    config=config,
                    colormap_section='viz_combined_cell_width',
                    title=f'Combined Great Lakes mesh: {name} outlines',
                    colorbar_label='cell width (km)',
                    plot_land=False,
                    extent=box,
                    patch_edge_color='black',
                    cell_indices=cell_indices,
                )
                span_km = _box_span_km(box)
                median_width_km = float(np.median(width_km[cell_indices]))
                zoom_lines.append(
                    f'{name}: nCells={cell_indices.size} '
                    f'span_km={span_km:.1f} '
                    f'median_width_km={median_width_km:.3f} '
                    f'span_km_over_median_width_km='
                    f'{span_km / median_width_km:.1f}'
                )

            basin_names = ds_topo.attrs['basin_names'].split(',')
            basin_index = ds_topo['basinIndex'].values
            basin_lines = [
                f'{name}: {int((basin_index == i).sum())}'
                for i, name in enumerate(basin_names)
            ]

            with open('debug_summary.txt', 'w', encoding='utf-8') as summary:
                summary.write(f'nCells: {ds_mesh.sizes["nCells"]}\n')
                summary.write(
                    f'cell_width_km_min: {float(width_km.min()):.6f}\n'
                )
                summary.write(
                    f'cell_width_km_median: {float(np.median(width_km)):.6f}\n'
                )
                summary.write(
                    f'cell_width_km_max: {float(width_km.max()):.6f}\n'
                )
                dc_edge_km = 1.0e-3 * ds_mesh['dcEdge'].values
                summary.write(
                    f'dc_edge_km_min: {float(dc_edge_km.min()):.6f}\n'
                )
                summary.write(
                    f'dc_edge_km_max: {float(dc_edge_km.max()):.6f}\n'
                )
                summary.write(
                    f'base_elevation_min: {float(base_elevation.min()):.6f}\n'
                )
                summary.write(
                    f'base_elevation_max: {float(base_elevation.max()):.6f}\n'
                )
                summary.write(
                    f'base_elevation_nan_count: '
                    f'{int(np.isnan(base_elevation.values).sum())}\n'
                )
                summary.write('basin_cell_counts:\n')
                for line in basin_lines:
                    summary.write(f'  {line}\n')
                summary.write(
                    'mesh_outline_zooms (span_km_over_median_width_km '
                    'should stay well under 250 for legible outlines):\n'
                )
                for line in zoom_lines:
                    summary.write(f'  {line}\n')


def _cell_width_km(ds_mesh):
    """
    Compute each cell's width in km as the mean of ``dcEdge`` over its
    edges.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        An MPAS mesh with ``dcEdge``, ``edgesOnCell``, and
        ``nEdgesOnCell``

    Returns
    -------
    width_km : numpy.ndarray
        The cell width in km, one value per cell
    """
    dc_edge_km = 1.0e-3 * ds_mesh['dcEdge'].values
    edges_on_cell = ds_mesh['edgesOnCell'].values
    n_edges_on_cell = ds_mesh['nEdgesOnCell'].values
    max_edges = edges_on_cell.shape[1]
    # edgesOnCell is one-based and zero-padded out to maxEdges, so the
    # padding is masked off rather than indexed
    valid = np.arange(max_edges)[None, :] < n_edges_on_cell[:, None]
    indices = np.where(valid, edges_on_cell - 1, 0)
    total = np.where(valid, dc_edge_km[indices], 0.0).sum(axis=1)
    return total / n_edges_on_cell


def _cells_in_box(ds_mesh, extent):
    """
    Find the indices of cells whose center falls within a lon/lat box.

    Parameters
    ----------
    ds_mesh : xarray.Dataset
        An MPAS mesh with ``lonCell`` and ``latCell`` in radians, with
        ``lonCell`` in ``[0, 2 * pi)``

    extent : tuple of float
        The (lon_min, lon_max, lat_min, lat_max) box, in degrees, with
        longitudes in ``[-180, 180)``

    Returns
    -------
    cell_indices : numpy.ndarray
        The indices of cells within the box
    """
    lon_min, lon_max, lat_min, lat_max = extent
    lon_deg = np.rad2deg(ds_mesh['lonCell'].values)
    lon_deg = np.where(lon_deg > 180.0, lon_deg - 360.0, lon_deg)
    lat_deg = np.rad2deg(ds_mesh['latCell'].values)
    mask = (
        (lon_deg >= lon_min)
        & (lon_deg <= lon_max)
        & (lat_deg >= lat_min)
        & (lat_deg <= lat_max)
    )
    return np.flatnonzero(mask)


def _pad_box(extent, pad):
    """Pad a (lon_min, lon_max, lat_min, lat_max) box by ``pad`` degrees"""
    lon_min, lon_max, lat_min, lat_max = extent
    return (lon_min - pad, lon_max + pad, lat_min - pad, lat_max + pad)


def _box_span_km(extent):
    """Approximate east-west span of a lon/lat box in km at its center"""
    lon_min, lon_max, lat_min, lat_max = extent
    mean_lat_rad = np.deg2rad(0.5 * (lat_min + lat_max))
    return (lon_max - lon_min) * 111.0 * np.cos(mean_lat_rad)
