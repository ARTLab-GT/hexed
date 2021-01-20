find_package(Catch2 REQUIRED)
add_executable(test)
target_sources(test PRIVATE
main.cpp
test_Basis.cpp
test_conservation.cpp
test_Grid.cpp
test_local.cpp
test_max_char_speed.cpp
test_neighbor.cpp
test_sinusoid.cpp
test_Solution.cpp
test_func.cpp
)
target_link_libraries(test PRIVATE cartdg tecio Catch2::Catch2)
