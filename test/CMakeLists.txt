find_package(Catch2 REQUIRED)
add_executable(unit)
target_sources(unit PRIVATE
# please maintain alphabetical order
testing_utils.cpp
test_Accessible_mesh.cpp
test_Basis.cpp
test_Boundary_condition.cpp
test_Characteristics.cpp
test_connection.cpp
test_Deformed_element.cpp
test_Derivative.cpp
test_erase_if.cpp
test_Element.cpp
test_Element_container.cpp
test_func.cpp
test_Face_permutation.cpp
test_filter_limit.cpp
test_Hanging_vertex_matcher.cpp
test_Interpreter.cpp
test_Iteration_status.cpp
test_Layer_sequence.cpp
test_math.cpp
test_Max_dt.cpp
test_Namespace.cpp
test_Occt.cpp
test_pde.cpp
test_Physical_basis.cpp
test_Prolong_refined.cpp
test_Restrict_refined.cpp
test_Row_index.cpp
test_Row_rw.cpp
test_read_csv.cpp
test_Simplex_geom.cpp
test_Solver.cpp
test_Stopwatch.cpp
test_Stopwatch_tree.cpp
test_Storage_params.cpp
test_Surface_geom.cpp
test_Tecplot_file.cpp
test_Tree.cpp
test_utils.cpp
test_Vertex.cpp
test_Vis_data.cpp
)
add_executable(interactive)
target_sources(interactive PRIVATE
test_otter_vis.cpp
test_Solver_vis.cpp
)
target_link_libraries(unit PRIVATE hexed Catch2::Catch2WithMain)
target_link_libraries(interactive PRIVATE hexed Catch2::Catch2WithMain)
