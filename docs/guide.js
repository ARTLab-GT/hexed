var guide =
[
    [ "Installation", "installation.html", [
      [ "Prerequisites", "installation.html#install_prereqs", [
        [ "OS", "installation.html#install_prereq_os", null ],
        [ "Python/pip", "installation.html#install_prereq_python", null ]
      ] ],
      [ "Install from Python package (the easy way)", "installation.html#install_py_pack", null ],
      [ "Building from source (the harder way)", "installation.html#install_source", [
        [ "Build options", "installation.html#build_options", [
          [ "source_dir", "installation.html#source_dir", null ],
          [ "build_dir", "installation.html#build_dir", null ],
          [ "n_build_procs", "installation.html#n_build_procs", null ],
          [ "venv", "installation.html#venv", null ],
          [ "compiler", "installation.html#compiler", null ],
          [ "verbose", "installation.html#verbose", null ],
          [ "use_system_paths", "installation.html#use_system_paths", null ],
          [ "use_env_paths", "installation.html#use_env_paths", null ],
          [ "internet", "installation.html#internet", null ],
          [ "build_mode", "installation.html#build_mode", null ],
          [ "max_row_size", "installation.html#max_row_size", null ],
          [ "threaded", "installation.html#threaded", null ],
          [ "n_threads", "installation.html#n_threads", null ],
          [ "profile", "installation.html#profile", null ],
          [ "use_xdmf", "installation.html#use_xdmf", null ],
          [ "use_tecio", "installation.html#use_tecio", null ],
          [ "build_tests", "installation.html#build_tests", null ],
          [ "run_tests", "installation.html#run_tests", null ],
          [ "test_args", "installation.html#test_args", null ],
          [ "gdb", "installation.html#gdb", null ],
          [ "build_docs", "installation.html#build_docs", null ],
          [ "obsessive_timing", "installation.html#obsessive_timing", null ],
          [ "build_wheel", "installation.html#build_wheel", null ],
          [ "install_wheel", "installation.html#install_wheel", null ],
          [ "architecture", "installation.html#architecture", null ],
          [ "sanitize", "installation.html#sanitize", null ]
        ] ],
        [ "Automatically-installed dependencies", "installation.html#auto_dependencies", null ]
      ] ]
    ] ],
    [ "Running Hexed", "running.html", null ],
    [ "Hexed Interface Language", "hil.html", [
      [ "TL;DR", "hil.html#hil_tldr", null ],
      [ "Introduction", "hil.html#Introduction", null ],
      [ "Types and Literals", "hil.html#literals", null ],
      [ "Operators", "hil.html#Operators", null ],
      [ "Variables and Assignment", "hil.html#variables", null ],
      [ "Macros", "hil.html#Macros", null ],
      [ "Builtin Variables", "hil.html#builtins", [
        [ "List of Builtin Variables", "hil.html#builtin_list", [
          [ "Physical/mathematical constants", "hil.html#hil_constants", null ],
          [ "huge", "hil.html#huge", null ],
          [ "nan", "hil.html#nan", null ],
          [ "newline", "hil.html#newline", null ],
          [ "false", "hil.html#false", null ],
          [ "true", "hil.html#true", null ],
          [ "loop", "hil.html#loop", null ],
          [ "print_emph", "hil.html#print_emph", null ],
          [ "print_type", "hil.html#print_type", null ],
          [ "repl", "hil.html#repl", null ],
          [ "quit", "hil.html#quit", null ],
          [ "ask", "hil.html#ask", null ],
          [ "except", "hil.html#except", null ]
        ] ],
        [ "exception", "hil.html#exception", [
          [ "exit", "hil.html#exit", null ],
          [ "throw", "hil.html#throw", null ],
          [ "system_time", "hil.html#system_time", null ],
          [ "steady_time", "hil.html#steady_time", null ]
        ] ]
      ] ],
      [ "Idioms", "hil.html#Idioms", [
        [ "Comments", "hil.html#Comments", null ],
        [ "Conditionals", "hil.html#Conditionals", null ],
        [ "Iteration", "hil.html#Iteration", null ],
        [ "Functions", "hil.html#Functions", null ]
      ] ]
    ] ],
    [ "Solver Parameters", "parameters.html", "parameters" ],
    [ "Boundary Conditions", "boundary_conditions.html", [
      [ "Freestream", "boundary_conditions.html#Freestream", null ],
      [ "Characteristic", "boundary_conditions.html#Characteristic", null ],
      [ "Pressure Outflow", "boundary_conditions.html#pressure_outflow", null ],
      [ "Outflow", "boundary_conditions.html#outflow", null ],
      [ "Nonpenetration", "boundary_conditions.html#Nonpenetration", null ],
      [ "No Slip", "boundary_conditions.html#no_slip", null ],
      [ "Expression", "boundary_conditions.html#expression", null ]
    ] ],
    [ "Error Types", "error_types.html", [
      [ "HIL Exception", "error_types.html#hil_exception", null ],
      [ "User Error", "error_types.html#user_error", null ],
      [ "Feature-Not-Implemented Error", "error_types.html#not_implemented_error", null ],
      [ "Internal Error", "error_types.html#internal_error", null ],
      [ "Numerical Exception", "error_types.html#numerical_error", null ],
      [ "Other Exceptions", "error_types.html#other_exceptions", null ],
      [ "Segmentation fault", "error_types.html#segfault", null ]
    ] ],
    [ "Mesh I/O", "mesh_io.html", null ],
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
        [ "Output", "geom_fitting.html#autotoc_md1", null ]
      ] ]
    ] ],
    [ "Contributing Guidelines", "contributing.html", [
      [ "General Procedures", "contributing.html#contrib_general", null ],
      [ "Pull Requests", "contributing.html#contrib_pull", null ],
      [ "Code Style", "contributing.html#contrib_style", [
        [ "File Formatting", "contributing.html#file_formatting", [
          [ "Brace Style", "contributing.html#brace_style", null ],
          [ "Class Style", "contributing.html#class_style", null ]
        ] ],
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
    ] ],
    [ "Weird Problems", "weird_problems.html", null ]
];