target_sources(hexed PRIVATE
# please maintain alphabetical order
Accessible_mesh.cpp
Basis.cpp
Boundary_condition.cpp
Boundary_func.cpp
Case.cpp
connection.cpp
criteria.cpp
Deformed_element.cpp
Domain_func.cpp
Element.cpp
Element_func.cpp
Equidistant.cpp
filter_limit.cpp
global_hacks.cpp
Hanging_vertex_matcher.cpp
History_monitor.cpp
hil_properties.cpp
Interpreter.cpp
iterative.cpp
Iteration_status.cpp
kernels_convective.cpp
kernels_diffusive.cpp
kernels_max_dt.cpp
Layer_sequence.cpp
Linearized.cpp
Linear_equation.cpp
Lock.cpp
math.cpp
Occt.cpp
Output_data.cpp
otter_vis.cpp
Physical_basis.cpp
Qpoint_func.cpp
read_csv.cpp
Simplex_geom.cpp
Solver.cpp
Solver_snap_faces.cpp
Spacetime_func.cpp
Stopwatch.cpp
Stopwatch_tree.cpp
Storage_params.cpp
Struct_expr.cpp
Surface_func.cpp
Surface_geom.cpp
stabilizing_art_visc.cpp
standard_atmosphere.cpp
Tecplot_file.cpp
Tree.cpp
thermo.cpp
Vertex.cpp
Vis_data.cpp
Xdmf_wrapper.cpp
)

add_executable(hexecute hexecute.cpp)
target_link_libraries(hexecute PRIVATE hexed)
add_executable(hil hil.cpp)
target_link_libraries(hil PRIVATE hexed)
