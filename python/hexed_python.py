import cppyy

cppyy.add_include_path("include")
cppyy.add_include_path("build_Release/include/auto")
cppyy.load_library("build_Release/libhexed")
cppyy.include("Solver.hpp")
status = cppyy.gbl.hexed.Iteration_status()
print(status.header())
#solver = cppyy.gbl.hexed.Solver(2, 6, 1.)
#print(solver.iteration_status().header())
