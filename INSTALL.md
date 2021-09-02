# Installation

## Quick Start
* Create and navigate to a build directory (git will ignore names that start with "build").
* Use `ccmake` to configure build options. This project has been configured to show you the available options and their defaults specifically with the `ccmake` interface for CMake.
  * Ensure that `TECIO_DIR` is set to a directory where `bin/libtecio.so` and `include/TECIO.h` can be found (the default *should* be correct for ARTLAB machines).
    Much to my disappointment, TecIO does not cooperate with CMake, so this option must be configured manually.
  * Set `CMAKE_INSTALL_PREFIX` to the location where you want to install CartDG. Alas, there is no reasonable default for this,
    because we do not have write permission above our respective home directories on ARTLAB machines.
  * Otherwise, the default options should be appropriate for a release build.
* `make install`
* CartDG is now ready to use. To verify the success of the build, execute the following:
  * `test/test`
    * This will execute the unit tests.
  * `python3 benchmark.py`
    * This will execute the performance benchmarking script.
  * `demo/demo`
    * This will run a very simple toy problem and output visualization files.

## Detailed instructions
This section is in progress. For now, please refer to Quick Start.
* Create install directory.
  * `cd`
  * `mkdir codes`
  * Edit bashrc and add:
    * `export PATH=$PATH:~/codes`
    * `export CXX_INCLUDE_PATH=$CXX_INCLUDE_PATH:~/codes/include`
* `pip3 install sympy`
* Install Catch2:
  * `git clone https://github.com/catchorg/Catch2.git`
  * `cd Catch2/`
  * `cmake -Bbuild -H. -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=~/codes`
  * `cmake --build build/ --target install`
* Install Eigen:
  * Download the Eigen [source code](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) (latest stable release).
  * Unpack (e.g. unzip) the Eigen source and copy it to `~/codes`. `~/codes` should now contain a directory named `eigen-X.X.X` depending
    on the current version. `cd` into this directory.
  * `cp -r Eigen/ ~/codes/include/`
    * `~/codes/include/` should have already been created when you installed Catch2, but if for whatever reason it does not exist, create it.
