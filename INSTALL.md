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
These instructions should work if you follow them exactly. If you are trying to adapt them for whatever reason, follow these general principles:
* If the instructions mention a directory that does not exist, create it.

Instructions:
* Create install directory.
  * `cd`
  * `mkdir codes`
  * Edit bashrc and add:
    * `export PATH=$PATH:~/codes`
    * `export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:~/codes/include`
* `pip3 install sympy`
* Install Catch2:
  * Visit the Catch2 [Releases](https://github.com/catchorg/Catch2/releases/tag/v2.13.7) page and download the source code
    for the latest one (do not clone the `devel` branch from the github repo, as that is not a stable release).
  * Unpack (e.g. unzip) the source code and place it in `~/codes`. `~/codes` should now contain a directory named something like `Catch2-2.xx.x`
    depending on what happens to be the current version. `cd` into this directory.
  * `cmake -Bbuild -H. -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=~/codes`
  * `cmake --build build/ --target install`
* Install Eigen:
  * Download the Eigen [source code](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) (latest stable release).
  * Unpack the Eigen source and copy it to `~/codes`. `~/codes` should now contain a directory named something like `eigen-X.X.X`.
    `cd` into this directory.
  * `cp -r Eigen/ ~/codes/include/`
  * build...
  * `python3 benchmark.py dim=2`
    * This runs the performance benchmarking script for 2 dimensions.
    * You should see three plots. One is titled "Fraction of time for 3-stage RK cycle: regular" and below
      that you should see something like "Total: 1.4e-06 s". Another is titled "Fraction... : deformed" and should say
      "Total: 1.2e-06 s". 
  * `python3 benchmark.py`
    * This runs the performance benchmarking script for 3 dimensions (the default).
    * You should see similar plots. The "regular" plot total should be 2.6e-5 s and the "deformed" should be 3.1e-5 s.
  * If any of these "total" numbers differ from this tutorial by more than a factor of 2 (better or worse) please
    report it to Micaiah. That said, you might try running if a few times first. Timing execution always involves some uncertainty,
    and sometimes there can be outlying results.
