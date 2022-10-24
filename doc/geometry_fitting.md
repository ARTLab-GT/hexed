# Geometry fitting
Hexed provides a script
[`hexed-fit-geom`](https://github.gatech.edu/ARTLab/hexed/blob/master/script/hexed-fit-geom)
for creating a smooth curve fit to a 2D geometry.
Given a geometry file, it computes parametric polynomial curve as a least squares fit to the geometry points.
A new geometry file is created that contains the smooth fit sampled at an arbitrary number of points, uniformly distributed in parameter space.
This is useful when a geometry file is available but the resolution and/or precision is insufficient.
The script is designed to work well for airfoils, but may also work for other geometries.

## Arguments
`hexed-fit-geom` accepts the following command line arguments in any order:
The following argument is required:
- **file name:** Name of the file containing the input data points to which a smooth curve will be fitted.
  The file may have a  `.geom` extension, in which case it must be in the [NASCART-GT](https://github.gatech.edu/ARTLab/NASCART-GT)
  geometry file format ([example](https://github.gatech.edu/ARTLab/NASCART-GT/blob/master/Examples/RAE2822_airfoil/RAE2822.geom)).
  Alternatively, it may have a `.txt` or `.csv` extension in which it must contain two columns of data representing *x* and  *y* coordinates
  delimited by a single space or comma, respectively.
  This geometry should be smooth.
  In particular, for airfoils, the coordinates should start and end at the trailing edge, not the leading.

The following arguments are optional:
- `--uncluster-trailing`: If this argument is provided, the distribution of data points (not output sample points) in parametric space will be modified
  to avoid overfitting near the trailing edge.
  Generally, for airfoils, geometry points should not be clustered near the trailing edge (although mesh points should).
  If they are, this option can help mitigate the resulting ill-conditioning of the curve fit.
  In any case, if you see the curve fit looking ugly near the trailing edge, consider trying this option.
- `--optimize`: If provided, distribution of data points in parametric space will be optimized to minimize the residual of the curve fit.
  This is designed to fix nonsmooth distributions of geometry points which can disrupt the polynomial interpolation.
  However, it is computationally expensive and not very robust, so only use it if necessary.
  If your curve fit is ugly and you notice rapid changes in the density of data points or irregular spacing, try this option.
- `--close`: If provided, the output geometry will be smoothly deformed such that the first and last points match exactly
  (since the curve fit will match them at best approximately).
  Otherwise, a segment will be added to connect the first and last points.
  For airfoils, it is generally recommended that you use this option, especially if the trailing edge of your original geometry was closed.
- `--quiet`: If provided, disable all command-line output and diagnostic plots.
- `--periodic`: Specifies that a trigonometric basis shall be used instead of a polynomial one.
  Use this if the end of your geometry curve is supposed to blend smoothly into the beginning (e.g., the starfish case).
  Do not use this if your geometry has an edge that coincides with the start/end of the curve (e.g., an airfoil).
- `--order=[n]`: Specifies an `[n]`th order curve fit.
  If not periodic, this is the number of polynomial coefficients.
  E.g. `--order=20` specifies a 19th-degree polynomial curve fit.
  If periodic, this is the number of cosine *and* sine functions that will be used.
  If not provided, defaults to 10.
- `--num-points=[n]`: Specifies that the output file should include `[n]` sample points.
  If not provided, defaults to 1000.

## Output
Depending on the `--quiet` argument, some diagnostic information may be printed and plots may be displayed.
An output geometry file will also be generated.

The first plot shows the curve fit in parametric space.
That is, the resulting curve is defined by the coordinates *(x(t), y(t))* where *t* is the abscissa of the plot and *x* and *y* are displayed on the ordinate.
This is the space in which the least-squares fit is being computed.
The distribution of sample points in parametric space can also be observed in this plot.
If you can get past the abstractness of it, it is useful for troubleshooting the curve fit.

The second plot shows the curve fit in physical space.
By this plot you can decide whether you're satisfied with the accuracy of the fit and the smoothness of the curve.

Finally, an output file named `[name]_fitted.geom` where `[name]` is the name of the input file with the extension removed will be generated.
This will contain the curve-fit geometry in the NASCART-GT 2D geometry file format.
