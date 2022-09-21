# Installation

## Quick Start
* Create and navigate to a build directory (git will ignore names that start with "build").
* Use `ccmake` to configure build options.
  This project has been configured to show you the available options and their defaults specifically with the `ccmake` interface for CMake.
  * The default options should be appropriate for a release build.
  * Note that `CMAKE_INSTALL_PREFIX` defaults to `~/.local` instead of `/usr/local`, so that it can be installed on the lab machines without root access.
* Ensure that your system can find Tecplot/TecIO header and library files. You may need to edit `CPLUS_INCLUDE_PATH`, `LIBRARY_PATH`, and `LD_LIBRARY_PATH`.
  Much to my disappointment, TecIO does not cooperate with CMake, so this must be setup manually.
* `make -j install`
* Hexed is now ready to use. To verify the success of the build, execute the following:
  * `test/unit` to run the unit tests.
  * This will execute the performance benchmarking script.
  * `demo/demo2d`. This will run a very simple toy problem and output tecplot visualization files.
    To create a prettier version of the visualization, execute `hexed-post-process ".*szplt"` and use Tecplot to view the file `all.lay`.

## Detailed instructions
This section is provided in case you need more information than the Quick Start section includes. It assumes you have access to an ARTLAB
machine and are familiar with basic directory and file management in Linux. If you have more than that, it is assumed that you can adapt these
instructions to suit your environment. If you have less than that, consider asking for help. If you have any trouble following these instructions,
please let Micaiah know so that this documentation can be improved.

### 1. Configure environment
When you compile and install programs that rely on other programs,
your system uses so-called *environment variables* to remember where all the different programs are kept.
It is common to install programs you compile in `/usr/local`.
However, this requires superuser privileges, which we do not 
* `cd`
* `mkdir codes`
* We need to tell the shell to look for libraries in `~/codes`. We do this by editing the file `~/.bashrc`. First,
  it is recommended that you make a backup copy:
  * `cp .bashrc .bashrc_backup`

  Now edit `.bashrc` and add the following lines at the end:
  ```bash
  TECPLOT_DIR=/opt/local/tecplot/360ex
  TECPLOT_LIB_PATH=$TECPLOT_DIR/bin/
  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$TECPLOT_DIR/include/
  export LIBRARY_PATH=$LIBRARY_PATH:$TECPLOT_LIB_PATH
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TECPLOT_LIB_PATH
  ```
    * If you are not on an ARTLAB machine, ensure that you have Tecplot installed and change this path to match your Tecplot installation.
* Close and reopen the terminal.

### 2. Install Python libraries
CartDG uses Python scripts for some small tasks during the build process and output post-processing and requires
a few 3rd-party libraries.
* Install PyTecplot (Tecplot's Python API):
  * `pip3 install pytecplot`
* Install sympy (used to generate quadrature rules):
  * `pip3 install sympy`
### 3. Install Catch2
* Visit the Catch2 [Releases](https://github.com/catchorg/Catch2/releases) page and download the source code
  for the latest one (do not clone the `devel` branch from the github repo, as that is not a stable release).
* Extract (e.g., unzip) the source code and place it in `~/codes`. `~/codes` should now contain a directory named something like `Catch2-2.XX.X`,
  except instead of "X"s there will be a version number.
* `cd ~/Catch2-2.XX.X` (replace the "X"s to match whatever the actual directory name is).
* `cmake -Bbuild -H. -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=~/codes`
* `cmake --build build/ --target install`
### 4. Install Eigen
* Download the Eigen [source code](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) (latest stable release).
* Unpack the Eigen source in `~/codes`. `~/codes` should now contain a directory named something like `eigen-X.X.X`, depending on the
  current version number.
* If `~/codes/include` does not exist, create it.
* `cd eigen-X.X.X` (replace "X"s with appropriate numbers).
* `cp -r Eigen/ ~/codes/include/`
### 5. Install CartDG
  * `cd`
    * If you would like to put CartDG somewhere other than `~` (e.g. `~/codes`, go there instead.
  * `git clone git@github.gatech.edu:ARTLab/CartDG.git`
  * `cd CartDG`
  * `mkdir build_Release`
  * `cd build_Release`
  * `ccmake ../`
    * This will open a sort-of-GUI for the tool CMake. This allows you to choose some options for how to build CartDG. You should see `EMPTY_CACHE`
      and no error messages.
    * Press `c`. You should now see a list of options.
    * Use the arrow keys to move the cursor and "Enter" to edit the option under the cursor. If you want to abort and exit the GUI, press `q`.
    * The only option you need to modify is `CMAKE_INSTALL_PREFIX`. Set that to `~/codes`. Editing other parameters is optional.
    * When you have set all the options, press `c`. Then press `g`. The GUI should exit.
  * `make -j install`
### 6. Verify success
*  You should still be in the build directory.
* `test/test`
  * This executes the unit tests. It should say "All tests passed".
* `demo/demo`
  * This runs a simple demonstration simulation. It should create a directory `simulation` containing a bunch of `.szplt` files.
* `cartdg-post-process simulation/`
  * This will post-process the `.szplt` files to create a file `simulation/all.lay`. Open this file with Tecplot, and you should see
    a vortex travel out one side of the domain and in the other. If you try to open the `.szplt` files directly, you should see the same
    thing but not as pretty.
* `python3 benchmark.py dim=2 n_side=200`
  * This runs the performance benchmarking script for 2 dimensions.
  * You should see three plots. One is titled "Fraction of time for 3-stage RK cycle: regular" and below
    that you should see something like "Total: 2.1e-06 s". Another is titled "Fraction... : deformed" and should say
    "Total: 3.0e-06 s". 
* `python3 benchmark.py`
  * This runs the performance benchmarking script for 3 dimensions (the default).
  * You should see similar plots. The "regular" plot total should be 1.8e-05 s and the "deformed" should be 3.1e-05 s.
* If any of these "total" numbers differ from this tutorial by more than a factor of 2 (better or worse) please
  report it to Micaiah. That said, you might try running if a few times first. Timing execution always involves some uncertainty,
  and sometimes there can be outlying results.
