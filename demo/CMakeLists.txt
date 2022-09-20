add_executable(demo2d demo2d.cpp)
add_executable(demo3d demo3d.cpp)
target_link_libraries(demo2d PRIVATE hexed)
target_link_libraries(demo3d PRIVATE hexed)
