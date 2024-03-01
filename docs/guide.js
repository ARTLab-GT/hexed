var guide =
[
    [ "Installation", "installation.html", [
      [ "Quick Start", "installation.html#quick", null ],
      [ "Detailed instructions", "installation.html#detailed", [
        [ "Step 1: Set up installation directory", "installation.html#configure", null ],
        [ "Step 2: install Python libraries", "installation.html#install_python", null ],
        [ "Step 3: install Eigen", "installation.html#install_eigen", null ],
        [ "Step 4: install HDF5", "installation.html#install_hdf5", null ],
        [ "Step 5: install Boost headers", "installation.html#install_boost", null ],
        [ "Step 6: install libxml2", "installation.html#install_xml", null ],
        [ "Step 7: install XDMF", "installation.html#install_xdmf", null ],
        [ "Step 8: install OpenCASCADE Technology", "installation.html#install_occt", null ],
        [ "Step 9: compile and install Hexed", "installation.html#install_hexed", null ],
        [ "Step 10: verify success", "installation.html#verify", null ],
        [ "Development build", "installation.html#dev_build", [
          [ "Step 1: Install additional python packages", "installation.html#additional_python", null ],
          [ "Step 2: install Catch2", "installation.html#install_catch2", null ],
          [ "Step 3: build unit tests", "installation.html#build_tests", null ],
          [ "Step 4: build in debug mode", "installation.html#build_debug", null ]
        ] ]
      ] ],
      [ "Troubleshooting", "installation.html#Troubleshooting", [
        [ "XDMF problems", "installation.html#xdmf_problems", null ]
      ] ]
    ] ],
    [ "Running Hexed", "running.html", null ],
    [ "Hexed Interface Language", "hil.html", [
      [ "Types and Literals", "hil.html#literals", null ],
      [ "Operators", "hil.html#Operators", null ],
      [ "Variables and Assignment", "hil.html#variables", null ],
      [ "Macros", "hil.html#Macros", null ],
      [ "Builtin Variables", "hil.html#builtins", [
        [ "List of Builtin Variables", "hil.html#builtin_list", [
          [ "Physical/mathematical constants", "hil.html#hil_constants", null ],
          [ "huge", "hil.html#huge", null ],
          [ "newline", "hil.html#newline", null ],
          [ "false", "hil.html#false", null ],
          [ "true", "hil.html#true", null ],
          [ "loop", "hil.html#loop", null ],
          [ "repl", "hil.html#repl", null ],
          [ "quit", "hil.html#quit", null ],
          [ "ask", "hil.html#ask", null ],
          [ "exit", "hil.html#exit", null ],
          [ "throw", "hil.html#throw", null ]
        ] ]
      ] ],
      [ "Idioms", "hil.html#Idioms", [
        [ "Comments", "hil.html#Comments", null ],
        [ "Conditionals", "hil.html#Conditionals", null ],
        [ "Iteration", "hil.html#Iteration", null ],
        [ "Functions", "hil.html#Functions", null ]
      ] ]
    ] ],
    [ "Solver Parameters", "parameters.html", [
      [ "Input Parameters", "parameters.html#inputs", [
        [ "Important Input Parameters", "parameters.html#important_inputs", [
          [ "n_dim", "parameters.html#n_dim", null ],
          [ "done", "parameters.html#done", null ],
          [ "local_time", "parameters.html#local_time", null ],
          [ "max_safety", "parameters.html#max_safety", null ],
          [ "max_time_step", "parameters.html#max_time_step", null ],
          [ "transport_model", "parameters.html#transport_model", null ],
          [ "thermal_bc", "parameters.html#thermal_bc", null ],
          [ "row_size", "parameters.html#row_size", null ],
          [ "attack", "parameters.html#attack", null ],
          [ "sideslip", "parameters.html#sideslip", null ],
          [ "art_visc_width", "parameters.html#art_visc_width", null ],
          [ "refine_if", "parameters.html#refine_if", null ],
          [ "mesh_extreme00, mesh_extreme01, mesh_extreme10, mesh_extreme11, mesh_extreme20, mesh_extreme21", "parameters.html#mesh_extreme", null ],
          [ "extremal_bc00, extremal_bc01, extremal_bc10, extremal_bc11, extremal_bc20, extremal_bc21", "parameters.html#extremal_bc", null ],
          [ "n_cheby_flow", "parameters.html#n_cheby_flow", null ],
          [ "n_cheby_av", "parameters.html#n_cheby_av", null ],
          [ "av_advect_iters", "parameters.html#av_advect_iters", null ],
          [ "flow_iters", "parameters.html#flow_iters", null ],
          [ "init_layer_splits", "parameters.html#init_layer_splits", null ],
          [ "max_ref_level", "parameters.html#max_ref_level", null ]
        ] ],
        [ "Useful Input Parameters", "parameters.html#useful_inputs", [
          [ "init_ref_level", "parameters.html#init_ref_level", null ],
          [ "surface_bc", "parameters.html#surface_bc", null ],
          [ "monitor_vars", "parameters.html#monitor_vars", null ],
          [ "temperature_offset", "parameters.html#temperature_offset", null ],
          [ "geometry_uncertainty", "parameters.html#geometry_uncertainty", null ],
          [ "art_visc_constant", "parameters.html#art_visc_constant", null ],
          [ "fix_therm_admis", "parameters.html#fix_therm_admis", null ],
          [ "working_dir", "parameters.html#working_dir", null ],
          [ "vis_freq", "parameters.html#vis_freq", null ],
          [ "write_freq", "parameters.html#write_freq", null ],
          [ "print_freq", "parameters.html#print_freq", null ],
          [ "unrefine_if", "parameters.html#unrefine_if", null ],
          [ "max_refinement_iters", "parameters.html#max_refinement_iters", null ],
          [ "n_smooth", "parameters.html#n_smooth", null ],
          [ "geom_n_segments", "parameters.html#geom_n_segments", null ],
          [ "vis_surface", "parameters.html#vis_surface", null ],
          [ "vis_field", "parameters.html#vis_field", null ],
          [ "vis_tecplot", "parameters.html#vis_tecplot", null ],
          [ "vis_xdmf", "parameters.html#vis_xdmf", null ],
          [ "vis_skew", "parameters.html#vis_skew", null ],
          [ "vis_field_vars", "parameters.html#vis_field_vars", null ],
          [ "vis_surface_vars", "parameters.html#vis_surface_vars", null ],
          [ "vis_n_sample", "parameters.html#vis_n_sample", null ],
          [ "print_vars", "parameters.html#print_vars", null ],
          [ "integrand_field", "parameters.html#integrand_field", null ],
          [ "integrand_surface", "parameters.html#integrand_surface", null ],
          [ "max_angle", "parameters.html#max_angle", null ],
          [ "max_deflection", "parameters.html#max_deflection", null ],
          [ "flood_fill_start0, flood_fill_start1, flood_fill_start2", "parameters.html#flood_fill_start", null ],
          [ "init_cond", "parameters.html#init_cond", null ],
          [ "fix_admis_max_safety", "parameters.html#fix_admis_max_safety", null ],
          [ "av_advect_max_safety", "parameters.html#av_advect_max_safety", null ],
          [ "av_diff_max_safety", "parameters.html#av_diff_max_safety", null ],
          [ "layer_split_points", "parameters.html#layer_split_points", null ],
          [ "diffusive_admissibility", "parameters.html#diffusive_admissibility", null ],
          [ "input_data", "parameters.html#input_data", null ]
        ] ],
        [ "Esoteric Input Parameters", "parameters.html#esoteric_inputs", [
          [ "mesh_extreme_eps", "parameters.html#mesh_extreme_eps", null ],
          [ "av_diff_iters", "parameters.html#av_diff_iters", null ],
          [ "av_diff_ratio", "parameters.html#av_diff_ratio", null ],
          [ "av_visc_mult", "parameters.html#av_visc_mult", null ],
          [ "av_unscaled_max", "parameters.html#av_unscaled_max", null ],
          [ "av_advect_max_res", "parameters.html#av_advect_max_res", null ],
          [ "vis_lts_constraints", "parameters.html#vis_lts_constraints", null ],
          [ "iter_width", "parameters.html#iter_width", null ],
          [ "heat_flux_coercion", "parameters.html#heat_flux_coercion", null ],
          [ "shock_sub_iters", "parameters.html#shock_sub_iters", null ],
          [ "ramping_initial_iters", "parameters.html#ramping_initial_iters", null ],
          [ "ramping_initial_safety", "parameters.html#ramping_initial_safety", null ],
          [ "ramping_final_safety", "parameters.html#ramping_final_safety", null ],
          [ "monitor_window", "parameters.html#monitor_window", null ],
          [ "monitor_samples", "parameters.html#monitor_samples", null ],
          [ "elementwise_art_visc", "parameters.html#elementwise_art_visc", null ],
          [ "elementwise_art_visc_pde", "parameters.html#elementwise_art_visc_pde", null ],
          [ "elementwise_art_visc_diff_ratio", "parameters.html#elementwise_art_visc_diff_ratio", null ]
        ] ]
      ] ],
      [ "Command Variables & Macros", "parameters.html#commands", [
        [ "List of Command Variables", "parameters.html#command_list", [
          [ "start", "parameters.html#start", null ],
          [ "mesh", "parameters.html#mesh", null ],
          [ "run", "parameters.html#run", null ],
          [ "iterate", "parameters.html#iterate", null ],
          [ "create_solver", "parameters.html#create_solver", null ],
          [ "init_refinement", "parameters.html#init_refinement", null ],
          [ "add_geom", "parameters.html#add_geom", null ],
          [ "refine", "parameters.html#refine", null ],
          [ "init_state", "parameters.html#init_state", null ],
          [ "visualize", "parameters.html#visualize", null ],
          [ "compute_residuals", "parameters.html#compute_residuals", null ],
          [ "compute_forces", "parameters.html#compute_forces", null ],
          [ "print_forces", "parameters.html#print_forces", null ],
          [ "plot_history", "parameters.html#plot_history", null ],
          [ "update", "parameters.html#update", null ],
          [ "integrate_field", "parameters.html#integrate_field", null ],
          [ "integrate_surface", "parameters.html#integrate_surface", null ],
          [ "split_layers", "parameters.html#split_layers", null ],
          [ "set_geometry_refinement", "parameters.html#set_geometry_refinement", null ],
          [ "set_shock", "parameters.html#set_shock", null ],
          [ "set_ramping", "parameters.html#set_ramping", null ],
          [ "set_aerothermo", "parameters.html#set_aerothermo", null ],
          [ "read_mesh", "parameters.html#read_mesh", null ],
          [ "read_state", "parameters.html#read_state", null ],
          [ "write_mesh", "parameters.html#write_mesh", null ],
          [ "write_state", "parameters.html#write_state", null ]
        ] ]
      ] ],
      [ "Output Parameters", "parameters.html#outputs", [
        [ "List of Output Parameters", "parameters.html#output_list", [
          [ "input_file", "parameters.html#input_file", null ],
          [ "version", "parameters.html#version", null ],
          [ "version_major", "parameters.html#version_major", null ],
          [ "version_minor", "parameters.html#version_minor", null ],
          [ "version_patch", "parameters.html#version_patch", null ],
          [ "commit", "parameters.html#commit", null ],
          [ "license", "parameters.html#license", null ],
          [ "iteration", "parameters.html#iteration", null ],
          [ "wall_time", "parameters.html#wall_time", null ],
          [ "flow_time", "parameters.html#flow_time", null ],
          [ "time_step", "parameters.html#time_step", null ],
          [ "header", "parameters.html#header", null ],
          [ "report", "parameters.html#report", null ],
          [ "residual_momentum", "parameters.html#residual_momentum", null ],
          [ "residual_density", "parameters.html#residual_density", null ],
          [ "residual_energy", "parameters.html#residual_energy", null ],
          [ "normalized_residual", "parameters.html#normalized_residual", null ],
          [ "art_visc_residual", "parameters.html#art_visc_residual", null ],
          [ "n_elements", "parameters.html#n_elements", null ],
          [ "performance_report", "parameters.html#performance_report", null ]
        ] ]
      ] ]
    ] ],
    [ "Boundary Conditions", "boundary_conditions.html", [
      [ "Freestream", "boundary_conditions.html#Freestream", null ],
      [ "Characteristic", "boundary_conditions.html#Characteristic", null ],
      [ "Pressure Outflow", "boundary_conditions.html#pressure_outflow", null ],
      [ "Outflow", "boundary_conditions.html#outflow", null ],
      [ "Nonpenetration", "boundary_conditions.html#Nonpenetration", null ],
      [ "No Slip", "boundary_conditions.html#no_slip", null ]
    ] ],
    [ "Error Types", "error_types.html", [
      [ "HIL Exceptions", "error_types.html#hil_exceptions", null ],
      [ "Numerical Errors", "error_types.html#numerical_errors", null ],
      [ "Other Exceptions", "error_types.html#other_exceptions", null ],
      [ "Segmentation faults", "error_types.html#segfaults", null ]
    ] ],
    [ "Mesh I/O", "mesh_io.html", [
      [ "Native file format", "mesh_io.html#mesh_format", [
        [ "Mesh representation", "mesh_io.html#mesh_repr", [
          [ "Elements", "mesh_io.html#Elements", null ],
          [ "Vertex deformation", "mesh_io.html#vertex_def", null ],
          [ "Face warping", "mesh_io.html#face_warping", null ],
          [ "Tree", "mesh_io.html#Tree", null ],
          [ "Connections", "mesh_io.html#Connections", null ]
        ] ],
        [ "File format", "mesh_io.html#file_format", null ]
      ] ],
      [ "Export formats", "mesh_io.html#mesh_export", null ]
    ] ],
    [ "Notation and conventions", "conventions.html", [
      [ "Units and physical quantities", "conventions.html#units", [
        [ "State vector", "conventions.html#state_vector", null ]
      ] ],
      [ "Terminology", "conventions.html#Terminology", null ],
      [ "Abbreviations", "conventions.html#Abbreviations", null ],
      [ "Storage order", "conventions.html#storage_order", null ]
    ] ],
    [ "Geometry fitting", "geom_fitting.html", [
      [ "Arguments", "geom_fitting.html#geom_fit_args", [
        [ "Output", "geom_fitting.html#autotoc_md0", null ]
      ] ]
    ] ],
    [ "Contributing guidelines", "contributing.html", [
      [ "General procedures", "contributing.html#contrib_general", null ],
      [ "Pull requests", "contributing.html#contrib_pull", null ],
      [ "Code style", "contributing.html#contrib_style", [
        [ "File formatting", "contributing.html#file_formatting", null ],
        [ "Names", "contributing.html#Names", null ],
        [ "Exceptions", "contributing.html#Exceptions", null ],
        [ "Miscellaneous", "contributing.html#Miscellaneous", null ]
      ] ]
    ] ],
    [ "Performance Benchmarking", "benchmarking.html", [
      [ "Overall performance", "benchmarking.html#overall", null ],
      [ "Kernel performance breakdown", "benchmarking.html#breakdown", [
        [ "inviscid NACA 0012 case", "benchmarking.html#naca0012", null ],
        [ "viscous flat plate case", "benchmarking.html#flat_plate", null ],
        [ "partial Blottner sphere case", "benchmarking.html#blottner_sphere", null ]
      ] ]
    ] ]
];