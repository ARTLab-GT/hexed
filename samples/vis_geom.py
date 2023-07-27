import hexed

#hexed.cpp.Occt_geom(hexed.cpp.Occt_geom.read("samples/wing_store.IGS"), 3).visualize("wing_store")
#hexed.cpp.Occt_geom(hexed.cpp.Occt_geom.read("samples/crmhl-2dcut.igs"), 2).visualize("multielement")
hexed.make_geom("samples/wing_store.STL").visualize("wing_store_stl")
hexed.cpp.Simplex_geom[3](hexed.cpp.Occt.triangles(hexed.cpp.Occt.read("samples/wing_store.IGS"))).visualize("wing_store_discretized")
