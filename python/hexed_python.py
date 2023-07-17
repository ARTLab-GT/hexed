import cppyy

path = "build_Debug"
cppyy.add_include_path("include")
cppyy.add_include_path(f"{path}/include/auto")
cppyy.load_library(f"{path}/libhexed")
cppyy.include("Solver.hpp")
#status = cppyy.gbl.hexed.Iteration_status()
#print(status.header())
solver = cppyy.gbl.hexed.Solver(2, 6, 1.)
#print(solver.iteration_status().header())
