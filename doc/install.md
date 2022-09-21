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
instructions according to your situation. If you have less than that, consider asking for help. If you have any trouble following these instructions,
please let Micaiah know so that this documentation can be improved.

### 1. Configure environment
When you compile and install programs that rely on other programs,
your system uses so-called *environment variables* to remember where all the different programs are kept.
It is common to install programs you compile in `/usr/local`.
However, this requires superuser privileges, which we do not have on the lab machines (I know, I know).
My solution is to install things in `~/.local` instead.
However, this requires us to edit an environment variable so that the system knows what we're doing.
Plus, we need to edit a *whole bunch* of environment variables to accommodate Tecplot being a little... *unique*.
* To start with, type `cd` to make sure you're in your home directory.
* We need to tell the shell to look for libraries in `~/codes`. We do this by editing the file `~/.bashrc`. First,
  it is recommended that you make a backup copy:
  * `cp .bashrc .bashrc_backup`
  * Now edit `.bashrc` and add the following lines at the end:
  ```bash
  TECPLOT_DIR=/opt/local/tecplot/360ex
  TECPLOT_LIB_PATH=$TECPLOT_DIR/bin/
  export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$TECPLOT_DIR/include/:~/.local/include/
  export LIBRARY_PATH=$LIBRARY_PATH:$TECPLOT_LIB_PATH
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TECPLOT_LIB_PATH
  ```
    * If you are not on an ARTLAB machine, ensure that you have Tecplot installed and change `TECIO_DIR` to match your Tecplot installation.
* Close and reopen the terminal.
* You also need somewhere to put all the code I'm about to make you download, as opposed to dumping it in your home directory.
  Thus, do `mkdir codes`.

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
* `cmake -Bbuild -H. -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=~/.local`
* `cmake --build build/ --target install`
### 4. Install Eigen
* Download the Eigen [source code](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) (latest stable release).
* Unpack the Eigen source in `~/codes`. `~/codes` should now contain a directory named something like `eigen-X.X.X`, depending on the
  current version number.
* `cd eigen-X.X.X` (replace "X"s with appropriate numbers).
* `cp -r Eigen/ ~/.local/include/`
### 5. Install Hexed
  * `cd ~/codes`
  * `git clone git@github.gatech.edu:ARTLab/hexed.git`
  * `cd hexed`
  * `mkdir build_Release`
  * `cd build_Release`
  * `ccmake ../`
    * This will open a sort-of-GUI for the tool CMake. This allows you to choose some options for how to build CartDG. You should see `EMPTY_CACHE`
      and no error messages.
    * Press `c`. You should now see a list of options. However, the default options are fine for a release build, so you don't have to do anything.
      * If you want to edit any of the options, use the arrow keys to move the cursor and Enter to edit the option under the cursor.
    * Press `c` again.
    * Press `g`. The GUI should exit.
  * `make -j install`
### 6. Verify success
*  You should still be in the build directory.
* `test/unit`
  * This executes the unit tests. It should say "All tests passed".
* `demo/demo2d`
  * This runs a simple demonstration simulation. It should create a bunch of `.szplt` files.
* `hexed-post-process "*.szplt"`
  * This will post-process the `.szplt` files to create a file `all.lay`. Open this file with Tecplot, and you should see
    a vortex travel out one side of the domain and in the other. If you try to open the `.szplt` files directly, you should see the same
    thing but not as pretty, because the quadrature points are overlayed on top of the visualization points and auxiliary variables have not been computed.
* If you want, you can also run `demo/demo3d`.
