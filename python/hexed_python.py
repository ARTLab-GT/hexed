import cppyy
from hexed_py_config import *

## \namespace hexed_python
# \brief Internal namespace for the Python API.
# \details In the Doxygen documentation,
# everything in the Python API will appear in the namespace `hexed_python` because it is defined in the file `python/hexed_python.py`
# However, when you actually use it, it lives in the module `hexed` (import it with `import hexed`).
# This distinction prevents naming ambiguities in the documentation,
# since both the C++ namespace and the Python module are called `hexed`.

## \namespace hexed_python.cpp
# \brief Namespace for accessing the C++ API in Python.
# \details This submodule provides Python bindings to the C++ namespace `hexed`
# via [cppyy](https://cppyy.readthedocs.io/en/latest/starting.html).
# For the most part, everything should work pretty much exactly as it would in C++.
# For example, to call `hexed::math::pow``(3, 2)` in Python you can do `hexed.cpp.math.pow(3, 2)`.
# However, please read the docs on cppyy if you want to do anything fancy.
# This submodule doesn't provide access to all the C++ definitions,
# but it should provide everything you need for high-level solver execution.

for path in include_paths:
    cppyy.add_include_path(path)
for lib in libraries:
    cppyy.load_library(lib)
cppyy.include("hexed/math.hpp")
cppyy.include("hexed/Solver_interface.hpp")
## \cond
cpp = cppyy.gbl.hexed
## \endcond
