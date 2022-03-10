find_package(Catch2 REQUIRED)
add_executable(test)
target_sources(test PRIVATE
# please maintain alphabetical order
main.cpp
testing_utils.cpp
test_Basis.cpp
test_conservation.cpp
test_deformed.cpp
test_Deformed_element.cpp
test_Deformed_grid.cpp
test_Derivative.cpp
test_Element.cpp
test_func.cpp
test_Ghost_boundary_condition.cpp
test_Hanging_vertex_matcher.cpp
test_Iteration_status.cpp
test_local.cpp
test_neighbor.cpp
test_math.cpp
test_observing.cpp
test_Physical_basis.cpp
test_Refined_face.cpp
test_Regular_grid.cpp
test_shock_capturing.cpp
test_sinusoid.cpp
test_Solution.cpp
test_Storage_params.cpp
test_Tecplot_file.cpp
test_Vertex.cpp
)
target_include_directories(test PRIVATE ${CMAKE_SOURCE_DIR}/kernels)
target_link_libraries(test PRIVATE cartdg kernels Catch2::Catch2)
