# Installation

## Quick Start
If look at these instructions and don't know what I'm talking about, skip to "Detailed instructions".
* Ensure that all of the following are installed:
  [Tecplot](https://www.tecplot.com/)/[TecIO](https://www.tecplot.com/products/tecio-library/),
  [Eigen](https://eigen.tuxfamily.org/), [Catch2](https://github.com/catchorg/Catch2),
  and the Python packages numpy, scipy, sympy, and pytecplot (for Python3).
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
* You can change environment variables by editing the file `.bashrc`.
  First, it is recommended that you make a backup copy:
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
* Install the numerical computing libraries numpy and scipy (althought you likely have them already):
  * `pip3 install numpy`
  * `pip3 install scipy`
* Install the symbolic computation library sympy (used to generate quadrature rules):
  * `pip3 install sympy`
* Install PyTecplot (Tecplot's Python API):
  * `pip3 install pytecplot`
### 3. Install Eigen
* Download the Eigen [source code](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) (latest stable release).
* Unpack the Eigen source in `~/codes`. `~/codes` should now contain a directory named something like `eigen-X.X.X`, depending on the
  current version number.
* `cd eigen-X.X.X` (replace "X"s with appropriate numbers).
* `cp -r Eigen/ ~/.local/include/`
### 4. Install Hexed
  * `cd ~/codes`
  * `git clone git@github.gatech.edu:ARTLab/hexed.git`
  * `cd hexed`
  * `mkdir build_Release`
  * `cd build_Release`
  * `ccmake ../`
    * This will open [ccmake](https://cmake.org/cmake/help/latest/manual/ccmake.1.html), which is a sort-of-GUI for [CMake](https://cmake.org/).
      This allows you to choose some options for how to build Hexed.
      You should see `EMPTY_CACHE` and no error messages.
    * Press `c`. You should now see a list of options. The default options are fine for a release build, so you don't have to do anything.
      * If you want to edit any of the options, use the arrow keys to move the cursor and Enter to edit the option under the cursor.
    * Press `c` again.
    * Press `g`. The GUI should exit (sometimes you have to hit `c` several times before it will let you exit with `g`).
      If something goes wrong and you can't get `g` to work, you can abort with `q`.
  * `make -j install`
### 5. Verify success
*  You should still be in the build directory.
* `demo/demo2d`
  * This runs a simple demonstration simulation. It should create a bunch of `.szplt` files.
* `hexed-post-process "*.szplt"`
  * This will post-process the `.szplt` files to create a file `all.lay`. Open this file with Tecplot, and you should see
    a vortex travel out one side of the domain and in the other. If you try to open the `.szplt` files directly, you should see the same
    thing but not as pretty, because the quadrature points are overlayed on top of the visualization points and auxiliary variables have not been computed.
* If you want, you can also run `demo/demo3d`.

## Development build
The above instructions are sufficient if you just want to run simulations.
However, if you're going to be doing any development of Hexed, I recommend you perform the following additional steps.

### 1. Install Catch2
Hexed uses Catch2 is for unit testing. Install Catch2 so that you can compile the unit tests.
* Visit the Catch2 [Releases](https://github.com/catchorg/Catch2/releases) page and download the source code
  for the latest one (do not clone the `devel` branch from the github repo, as that is not a stable release).
* Extract (e.g., unzip) the source code and place it in `~/codes`. `~/codes` should now contain a directory named something like `Catch2-2.XX.X`,
  except instead of "X"s there will be a version number.
* `cd ~/codes/Catch2-2.XX.X` (replace the "X"s to match whatever the actual directory name is).
* `cmake -Bbuild -H. -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=~/.local`
* `cmake --build build/ --target install`

### 2. Build unit tests
Now you need to update your Hexed build to include the unit tests, as follows:
* `cd ~/codes/hexed/build_Release`
* `ccmake ..`
* Set the option `build_tests` to `ON`.
  You can accomplish this by using the arrow keys to move the cursor down to the appropriate line and then hitting "Enter" to toggle between `OFF` and `ON`.
* Hit `c` as many times as you need to and then `g`.
* `make -j` (you only need `install` again if you've changed any of the code and want to propogate the new version to NASCART-GT).
* The code is now compiled. you can run the unit tests with the command `test/unit` (executed from the build directory, still).
* You should see "All tests passed".

### 3. Create debug build
So far, we've only compiled Hexed in Release mode, which is designed to run as fast as possible at the expense of ease of debugging.
That's what you want when you're running simulations, but if you're in the process of implementing and debugging new features,
you would probably prefer to be able to use debugging tools like [GDB](https://www.sourceware.org/gdb/) and the [sanitizers](https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html).
So, we will create an new build directory where we compile in Debug mode.
* `cd ~/codes/hexed`
* `mkdir build_Debug`
* `ccmake ..`.
  * this time, set `CMAKE_BUILD_TYPE` to `Debug` ("Enter" toggles between `Release` and `Debug`),
    `build_tests` to `ON`, and `sanitize` to `ON`.
  * `c` as many times as you need and then `g`.
* `make -j` (**don't** add `install`, since you still want to be running simulations with your nice and fast Release mode).
* `test/unit`. Again, you should see "All tests passed".
   However, now if you were to, for example, write to an array out of bounds, you will get an error message with the line number of the problem
   instead of just a segfault.
