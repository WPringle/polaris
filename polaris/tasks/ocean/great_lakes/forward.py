from polaris.ocean.model import OceanModelStep, get_time_interval_string


class Forward(OceanModelStep):
    """
    A step for performing a short forward smoke-test run on the Great
    Lakes combined mesh and idealized initial condition.

    This validates that the combined mesh, vertical grid, and idealized
    initial condition actually run in Omega/MPAS-Ocean -- it is not a
    tuned, physically realistic simulation (see
    docs/design_docs/great_lakes_mesh.md).

    Attributes
    ----------
    init : polaris.tasks.ocean.great_lakes.init.Init
        The step that builds the mesh/vert_coord/init files this step
        runs from
    """

    def __init__(self, component, subdir, config, init):
        """
        Create the step

        Parameters
        ----------
        component : polaris.tasks.ocean.Ocean
            The ocean component the step belongs to

        subdir : str
            The subdirectory for the step

        config : polaris.config.PolarisConfigParser
            A shared config parser

        init : polaris.tasks.ocean.great_lakes.init.Init
            The step that builds the mesh/vert_coord/init files
        """
        super().__init__(
            component=component,
            name='forward',
            subdir=subdir,
            openmp_threads=1,
            graph_target=f'{init.path}/combined_graph.info',
            update_eos=True,
        )
        self.set_shared_config(config, link='great_lakes_forward.cfg')
        self.init = init

        self.add_horiz_mesh_input_file(work_dir_target=f'{init.path}/mesh.nc')
        self.add_vert_coord_input_file(
            work_dir_target=f'{init.path}/vert_coord.nc'
        )
        self.add_init_input_file(work_dir_target=f'{init.path}/init.nc')

        self.add_yaml_file('polaris.ocean.config', 'output.yaml')
        self.add_output_file(filename='output.nc')

        for baseline, time_index_end in [('init', -1)]:
            self.add_property_check(
                filename='output.nc',
                check_properties=[
                    'mass conservation',
                    'salt conservation',
                    'energy conservation',
                ],
                baseline=baseline,
                time_index_end=time_index_end,
            )

    def compute_cell_count(self):
        """
        The number of cells in the combined mesh, used to constrain
        resources. The combined mesh is fixed (not a family of
        resolutions), so this is a config option rather than a formula.

        Returns
        -------
        cell_count : int
            The approximate number of cells in the mesh
        """
        return self.config.getint('great_lakes_forward', 'approx_cell_count')

    def dynamic_model_config(self, at_setup):
        """
        Add model config options and yaml files that depend on config
        options, both at setup and runtime

        Parameters
        ----------
        at_setup : bool
            Whether this method is being run during setup of the step, as
            opposed to at runtime
        """
        super().dynamic_model_config(at_setup)

        config = self.config
        section = config['great_lakes_forward']

        btr_dt = section.getfloat('btr_dt')
        dt_btr_ratio = section.getfloat('dt_btr_ratio')
        dt = btr_dt * dt_btr_ratio
        del2 = section.getfloat('del2')
        del4 = section.getfloat('del4')
        disable_del4 = section.getboolean('disable_del4')
        bottom_drag_coeff = section.getfloat('bottom_drag_coeff')
        run_duration_hours = section.getfloat('run_duration_hours')
        output_interval_hours = section.getfloat('output_interval_hours')

        model = config.get('ocean', 'model')
        time_integrator = (
            'SplitExplicitRK2' if model == 'omega' else 'split_explicit'
        )

        replacements = dict(
            dt=get_time_interval_string(seconds=dt),
            btr_dt=get_time_interval_string(seconds=btr_dt),
            run_duration=get_time_interval_string(
                seconds=run_duration_hours * 3600.0
            ),
            time_integrator=time_integrator,
            del2=f'{del2:g}',
            del4=f'{del4:g}',
            use_del4=str(not disable_del4).lower(),
            bottom_drag_coeff=f'{bottom_drag_coeff:g}',
            output_interval=get_time_interval_string(
                seconds=output_interval_hours * 3600.0
            ),
            output_freq=f'{output_interval_hours:g}',
            output_freq_units='hours',
        )

        self.add_yaml_file(
            'polaris.tasks.ocean.great_lakes',
            'forward.yaml',
            template_replacements=replacements,
        )
