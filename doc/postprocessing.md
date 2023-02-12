# Visualization
Hexed can output Tecplot binary files (`.szplt`). These can be viewed in two ways.
You can simply open them normally in Tecplot, or you can use the postprocessing script (see below) to generate a layout file.
If you choose to open them normally, you will see the following:
- The domain will be divided into a zone for every mesh element.
- Each element/zone will contain a mesh of uniformly spaced sample points.
  These sample points are chosen to provide a clear view of the solution polynomial and do not represent the actual degrees of freedom.
- The dataset will contain the variables [`state0`, `state1`, ...].
  These are the conserved variables in the order [momentum, density, total energy].
  For example, in 2D, `state0` is x-momentum, `state1` is y-momentum, `state2` is density, and `state3` is total energy.
  Depending on whether you used artificial viscosity, the variable `artificial_viscosity_coefficient` may also be present.

Use the "Edge" feature to view the element mesh.
You are advised not to use the "Mesh" feature because it will show you the sample points, which you probably don't care about.
To visualize flow variables, use the "Calculate Variables..." function in the "Analyze" tab to compute quantities of practical interest (e.g. pressure, Mach number).
Since doing that every time you open a file would be a pain, Hexed also provides a Tecplot script to compute some of the popular ones for you.
To run it, select "Play Macro/Script..." in the "Scripting" tab, which should open a file browser.
Navigate to the Hexed build directory and select `format.mcr`.
You should immediately see a contour plot of Mach number.

## Postprocessing script
The [post-processing script](https://github.gatech.edu/ARTLab/hexed/blob/master/script/hexed-wrapped-post-process.in)
can be executed with the command `hexed-post-process`.
Assuming you [properly installed](install.md) Hexed,
you should be able to use this command from any directory.

### Arguments
The script takes two arguments, both optional, in arbitrary order:
- A [regular expression](https://docs.python.org/3/library/re.html#regular-expression-syntax) specifying the files to process.
  If you only want to process one file, you can just use the name of the file.
  If you want to view multiple files at once (for example, to see the transient history of a simulation),
  use python regex syntax (not a shell wildcard!) and enclose it in quotes so that the shell won't get upset about the special characters (examples below).
  All files specified must be Tecplot subzone data files (`.szplt`).
  If this argument is not specified, it defaults to `.*\.szplt` (all files in your working directory with a `.szplt` extension).
- `colormap=` followed by the name of the colormap you want to for the contour plots.
  Default: `colormap=hexed_plasma`.
  Of course, you can change this once you open Tecplot -- the script only affects what the colormap will be by default.
  Before setting the colormap, the script will import some of the
  [matplotlib colormaps](https://matplotlib.org/stable/gallery/color/colormap_reference.html)
  (namely as many of the perceptually uniform sequential colormaps as are available).
  These can be specified as the name of the colormap in matplotlib prefixed with `hexed_`.
  The ones available on the ARTLab machines are
  - inferno
  - magma
  - plasma
  - viridis

  All of these colormaps will be loaded, but when you open Tecplot it may or may not remember the ones you didn't choose.
  You can also import these colormaps yourself from the `.map` files in the `colormaps` subdirectory of you build directory.
  Why go through all this trouble when Tecplot already has plenty of colormaps to choose from?
  Because I personally hate every single one of the default Tecplot colormaps.
  As this [page](https://matplotlib.org/stable/tutorials/colors/colormaps.html) from the matplotlib documentation explains,
  an intuitive mapping (and viable conversion to grayscale) is achieved by associating the lightness of the color with the value of the data.
  All Tecplot provides is some rainbow-style colormaps which spectacularly fail to do this,
  and some sequential colormaps that have just about zero contrast.
  They even provide a colormap called "magma", but unlike the usual magma colormap both ends of the scale are light and the middle is dark.
  It's as though they went out of their way to deprive the user of legible colormaps.
  That said, if you want to specify one of these colormaps *anyway*, you can.
  Since most of the default colormap names have spaces, two consecutive underscores in the colormap argument will be converted to a space (see examples).
  
  tl;dr: use `colormap=hexed_[colormap name]` to use one of the matplotlib colormaps.
  
### Output
The script outputs a Tecplot layout file `all.lay`.
When you're ready to look at your data, simply open this file in Tecplot.
When interpreting this data, note:
- "mesh" may or may not show you anything useful
- "contour" shows you the solution polynomials sampled a mesh of uniformly spaced (in the reference element) points
  (usually there will be more sample points than quadrature points).
  In my view, this way of visualizing the solution most accurately reflects the math of the DG method.
  However, don't confuse the sample points with the quadrature points where the solution data is actually stored.
- "edge" shows you the edges of the elements.
  This is the easiest way to look at the mesh, but if you're used to looking at finite volume meshes,
  remember that each element contains **much** more information than a finite volume cell, so the meaning of this mesh is somewhat different.
  
### Examples
- `hexed-post-process`
  - visualize all `.szplt` files in your working directory with the matplotlib colormap "plasma".
- `hexed-post-process computational_surface_000052200.szplt`
  - visualize the file `computational_surface_000052200.szplt` with the plasma colormap.
- `hexed-post-process "flowfield.*" colormap=hexed_viridis`
  - visualize all files with names that start with `flowfield` using the matplotlib colormap "viridis".
- `hexed-post-process colormap=Large__Rainbow`
  - visualize all `.szplt` files in your workign directory with the Tecplot builtin colormap "Large Rainbow".
