# Installation

## Quick Start
* Create and navigate to a build directory (git will ignore names that start with "build").
* Use `ccmake` to configure build options. This project has been configured to show you the available options and their defaults specifically with the `ccmake` interface for CMake.
  * Ensure that `TECIO_DIR` is set to a directory where `bin/libtecio.so` and `include/TECIO.h` can be found (the default *should* be correct for ARTLAB machines).
    Much to my disappointment, TecIO does not cooperate with CMake, so this option must be configured manually.
  * Set `CMAKE_INSTALL_PREFIX` to the location where you want to install CartDG. Alas, there is no reasonable default for this,
    because we do not have write permission above our respective home directories on ARTLAB machines.
  * Otherwise, the default options should be appropriate for a release build.
* `make -j install`
* CartDG is now ready to use. To verify the success of the build, execute the following:
  * `test/test`
    * This will execute the unit tests.
  * `python3 benchmark.py`
    * This will execute the performance benchmarking script.
  * `demo/demo`
    * This will run a very simple toy problem and output visualization files.

## Detailed instructions
This section is provided in case you need more information than the Quick Start section includes. It assumes you have access to an ARTLAB
machine and are familiar with managing directories and files in Linux. If you have more than that, it is assumed that you can adapt these
instructions to suit your environment. If you have less than that, consider asking for help. If you have any trouble following these instructions,
please discuss it so that this documentation can be improved.
1. Create install directory. You need somewhere to install libraries so that your other codes can find them. On a system where you have
   superuser privileges, you can use `/usr/local`. However, on the ARTLAB machines, we can only write within our home directory, so we need
   to create our own install directory. If you already have one set up, you can use that. Otherwise, follow these instructions to create one:
  * `cd`
  * `mkdir codes`
  * We need to tell the shell to look for libraries in `~/codes`. Edit `~/.bashrc` and add the following lines at the end:
    * `export PATH=$PATH:~/codes/bin`
    * `TECPLOT_DIR=/opt/local/tecplot/360ex`
      * If you are not on an ARTLAB machine, change this path to match your Tecplot installation.
    * `TECPLOT_LIB_PATH=$TECPLOT_DIR/bin/`
    * `export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:~/codes/include/:$TECPLOT_DIR/include/`
    * `export LIBRARY_PATH=$LIBRARY_PATH:$TECPLOT_LIB_PATH`
    * `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TECPLOT_LIB_PATH`
    * Close and reopen the terminal.
2. Install PyTecplot (Tecplot's Python API):
  * `pip3 install pytecplot`
2. Install sympy (used by Cartdg to generate quadrature rules):
  * `pip3 install sympy`
3. Install Catch2:
  * Visit the Catch2 [Releases](https://github.com/catchorg/Catch2/releases) page and download the source code
    for the latest one (do not clone the `devel` branch from the github repo, as that is not a stable release).
  * Unpack (e.g. unzip) the source code and place it in `~/codes`. `~/codes` should now contain a directory named something like `Catch2-2.XX.X`
    depending on what happens to be the current version. `cd` into this directory.
  * `cmake -Bbuild -H. -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=~/codes`
  * `cmake --build build/ --target install`
4. Install Eigen:
  * Download the Eigen [source code](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) (latest stable release).
  * Unpack the Eigen source in `~/codes`. `~/codes` should now contain a directory named something like `eigen-X.X.X`.
    `cd` into this directory.
  * `cp -r Eigen/ ~/codes/include/`
    * If `~/codes/include` does not exist, create it.
5. Install CartDG:
  * Go to whatever directory you want to put CartDG in (`~` or `~/codes` works).
  * `git clone git@github.gatech.edu:ARTLab/CartDG.git`
  * `cd CartDG`
  * `mkdir build_Release`
  * `cd build_Release`
  * `ccmake ../`
    * This will open a GUI for the tool CMake. This allows you to choose some options for how to build CartDG. You should not see any errors.
    * Press `c`. You should now see a list of options.
    * Use the arrow keys to move the cursor and "Enter" to edit the option under the cursor.
    * If you are on an ARTLAB machine, the only option you need to edit is `CMAKE_INSTALL_PREFIX`. Set that to `~/codes`.
    * If you are not on an ARTLAB machine, you will need to also edit the `TECIO_DIR` option (and install Tecplot if you don't already have it).
    * Leave the other options. If you are developing CartDG, consider making another directory `build_Debug` where you set `CMAKE_BUILD_TYPE` to
      `Debug` and `sanitize` to `ON`.
    * When you have set all the options, press `c`. Then press `g`. The GUI should exit.
  * `make -j install`
6. Verify that everything works:
  * `test/test`
    * This executes the unit tests. It should say "All tests passed".
  * `demo/demo`
    * This runs a simple demonstration simulation. It should generate a bunch of `.szplt` files. If you open them in Tecplot, you should see
      a vortex travel out one side of the domain and come back on the other.
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
