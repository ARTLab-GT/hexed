find_package(Catch2 REQUIRED)
add_executable(unit)
target_sources(unit PRIVATE
# please maintain alphabetical order
main.cpp
testing_utils.cpp
test_Accessible_mesh.cpp
test_Basis.cpp
test_Boundary_condition.cpp
test_connection.cpp
test_Deformed_element.cpp
test_Derivative.cpp
test_erase_if.cpp
test_Element.cpp
test_Element_container.cpp
test_func.cpp
test_Face_permutation.cpp
test_hll.cpp
test_Hanging_vertex_matcher.cpp
test_Iteration_status.cpp
test_Local_cartesian.cpp
test_Local_deformed.cpp
test_math.cpp
test_Mcs_cartesian.cpp
test_Mcs_deformed.cpp
test_Physical_basis.cpp
test_Prolong_refined.cpp
test_Refined_face.cpp
test_Restrict_refined.cpp
test_Solver.cpp
test_Stopwatch.cpp
test_Stopwatch_tree.cpp
test_Storage_params.cpp
test_Surface_rotation.cpp
test_Tecplot_file.cpp
test_Vertex.cpp
test_Vis_data.cpp
test_Write_face.cpp
)
target_link_libraries(unit PRIVATE cartdg Catch2::Catch2)
add_executable(interactive)
target_sources(interactive PRIVATE
main.cpp
test_otter_vis.cpp
test_Solver_vis.cpp
)
target_link_libraries(interactive PRIVATE cartdg Catch2::Catch2)
