#ifndef HEXED_SOLVER_HPP_
#define HEXED_SOLVER_HPP_

#include "Accessible_mesh.hpp"
#include "Gauss_legendre.hpp"
#include "Spacetime_func.hpp"
#include "Element_func.hpp"
#include "Iteration_status.hpp"
#include "Stopwatch_tree.hpp"
#include "config.hpp"
#include "kernel_factory.hpp"
#include "Transport_model.hpp"
#if HEXED_USE_OTTER
#include <otter/plot.hpp>
#include <otter/colormap.hpp>
#include <otter/colors.hpp>
#endif

namespace hexed
{

/*! \brief The main class that basically runs everything.
 * \details If you want to run a simulation with hexed, your workflow should be roughly the following:
 * -# construct a `Solver` object
 * -# interact with the `mesh()` object to build the mesh topology
 * -# call `snap_vertices()` and then `relax_vertices()` repeatedly until you have an acceptable body-fitted mesh.
 * -# call `snap_faces()` to get a good high-order surface fit
 * -# call `calc_jacobian()` to initialize internal parameters based on the mesh
 * -# call `initialize()` to initialize the flow state
 * -# call `update()` repeatedly to progress the simulation
 * -# call the some of the functions in output section to get the data you want from the simulation
 */
class Solver
{
  Storage_params params;
  Accessible_mesh acc_mesh;
  Gauss_legendre basis;
  Iteration_status status;
  Stopwatch_tree stopwatch;
  bool use_art_visc;
  bool fix_admis;
  int av_rs;
  std::unique_ptr<Kernel<Element&>> write_face;
  bool is_local_time;
  Transport_model visc;
  Transport_model therm_cond;
  int last_fix_vis_iter = std::numeric_limits<int>::min();

  void share_vertex_data(Element::vertex_value_access, Vertex::reduction = Vertex::vector_max);
  void fix_admissibility(double stability_ratio);
  void apply_state_bcs();
  void apply_flux_bcs();
  void apply_avc_diff_bcs();
  void apply_avc_diff_flux_bcs();
  void apply_fta_flux_bcs();
  void compute_inviscid(double dt, int i_stage);
  void compute_viscous(double dt, int i_stage);
  void compute_fta(double dt, int i_stage);
  void compute_advection(double dt, int i_stage);
  void compute_avc_diff(double dt, int i_stage);
  void fta(double dt, int i_stage);
  bool use_ldg();
  double max_dt();

  public:
  /*! \name scheme parameters
   * Feel free to tweak these at runtime to influence the behavior of the solver
   */
  //!\{
  double fix_admis_stab_rat = .7; //!< staility ratio for fixing thermodynamic admissibility.
  //!\}
  /*! \name artificial viscosity parameters
   * Tweakable parameters specifically affecting smoothness-based artificial viscosity calculation.
   * You are permitted to mess with these dynamically at runtime.
   */
  //!\{
  double av_diff_ratio = 5e-3; //!< ratio of diffusion time to advection width
  double av_visc_mult = 30.; //!< final scaling parameter applied to artificial viscosity coefficient
  double av_unscaled_max = 2e-3; //!< maximum artificial viscosity coefficient before scaling (i.e. nondimensional)
  double av_advect_stab_rat = .2; //!< stability ratio for advection
  double av_diff_stab_rat = .5; //!< stability ratio for diffusion
  int av_advect_iters = 2; //!< number of advection iterations to run each time `set_art_visc_smoothness` is called
  int av_diff_iters = 1; //!< number of diffusion iterations to run each time `set_art_visc_smoothness` is called
  //!\}

  /*!
   * \param n_dim number of dimensions
   * \param row_size row size of the basis (see \ref Terminology)
   * \param root_mesh_size sets the value of `Mesh::root_mesh_size()`
   * \param local_time_stepping whether to use local or global time stepping
   * \param viscosity_model determines whether the flow has viscosity (natural, not artificial) and if so, how it depends on temperature
   * \param thermal_conductivity_model determines whether the flow has thermal conductivity and if so, how it depends on temperature
   * \details If `viscosity_model` and `thermal_conductivity_model` are both `inviscid` _and_ you don't turn on artificial viscosity,
   * you will be solving the pure inviscid flow equations.
   * Otherwise, you will be solving the viscous flow equations using the LDG scheme,
   * potentially with some of the diffusion coefficients
   * (artificial viscosity, natural viscosity, thermal conductivity) set to zero.
   */
  Solver(int n_dim, int row_size, double root_mesh_size, bool local_time_stepping = false,
         Transport_model viscosity_model = inviscid, Transport_model thermal_conductivity_model = inviscid);
  virtual ~Solver() = default;

  //! \name setup
  //!\{
  /*! \brief fetch the `Mesh`.
   * An object the user can use to build the mesh.
   * Note that whenever elements are added, the flow state, and Jacobian are uninitialized,
   * the time step scale is uniformly 1,
   * and the mesh quality may be poor.
   * The functions below must be used to complete the setup
   * before any flow calculation can begin.
   */
  inline Mesh& mesh() {return acc_mesh;}
  inline Storage_params storage_params() {return params;}
  //! moves all vertices to the mean of the current position and the mean of the neighbors' positions
  void relax_vertices(double factor = .5);
  /*! \brief apply `Mesh_bc`s
   * \note if some elements participate in multiple BCs, then snapping may not satisfy all exactly.
   *       However, if performed multiple times, it should converge;
   */
  void snap_vertices();
  //! warps the boundary elements such that the element faces coincide with the boundary at their quadrature points.
  void snap_faces();
  /*!
   * compute the Jacobian of all elements based on the current position of the vertices
   * and value of any face warping.
   * Mesh topology must be valid (no duplicate or missing connections) before calling this function.
   */
  void calc_jacobian();
  //! set the flow state
  void initialize(const Spacetime_func&);
  void set_art_visc_off(); //!< turns off artificial viscosity
  void set_art_visc_constant(double); //!< turns on artificial viscosity and initializes coefficient to a uniform value
  void set_art_visc_row_size(int); //!< modify the polynomial order of smoothness-based artificial viscosity (must be <= row size of discretization (which is the default))
  void set_fix_admissibility(bool); //!< turns on/off the thermodynamic admissibility-preserving scheme (increases robustness at some computational overhead)
  /*!
   * set `Element::resolution_badness` for each element according to `func`.
   * Resoltuion badness can be evaluated via `sample(ref_level, is_deformed, serial_n, Resolution_badness())`.
   * This function does some additional work to enforce some conditions on the resolution badness of neighboring elements.
   * Thus, use this function rather than just `sample(ref_level, is_deformed, serial_n, func)` directly.
   */
  void set_resolution_badness(const Element_func& func);
  //! set resolution badness of each element to be at least the maximum badness of any elements extruded from it
  void synch_extruded_res_bad();
  //!\}

  //! \name time marching
  //!\{
  /*!
   * March the simulation forward by a time step equal to `stability_ratio`
   * times the estimated maximum stable time step.
   */
  void update(double stability_ratio = 0.8);
  bool is_admissible(); //!< check whether flowfield is admissible (e.g. density and energy are positive)
  void set_art_visc_smoothness(double advect_length); //!< updates the aritificial viscosity coefficient based on smoothness of the flow variables
  //! an object providing all available information about the status of the time marching iteration.
  Iteration_status iteration_status();
  /*!
   * reset any variables in `iteration_status()` that count something since the last call to `reset_counters()`
   * (e.g. number of artificial viscosity iterations)
   */
  void reset_counters();
  //!\}

  //! \name output
  //! functions that compute some form of output data
  //!\{
  std::vector<double> sample(int ref_level, bool is_deformed, int serial_n, int i_qpoint, const Qpoint_func&); //!< evaluate arbitrary functions at arbitrary locations
  std::vector<double> sample(int ref_level, bool is_deformed, int serial_n, const Element_func&);
  //! obtain performance data
  inline const Stopwatch_tree& stopwatch_tree() {return stopwatch;}
  //! compute an integral over the entire flow field at the current time
  std::vector<double> integral_field(const Qpoint_func& integrand);
  //! compute an integral over all surfaces where a particular boundary condition has been enforced
  std::vector<double> integral_surface(const Boundary_func& integrand, int bc_sn);
  /*! compute the min and max of variables over entire flow field. layout: `{{var0_min, var0_max}, {var1_min, var1_max}, ...}`
   * bounds are approximated by uniformly sampling a block `n_sample`-on-a-side in each element
   */
  std::vector<std::array<double, 2>> bounds_field(const Qpoint_func&, int n_sample = 20);
  #if HEXED_USE_TECPLOT
  //! write a visualization file describing the entire flow field (but not identifying surfaces)
  void visualize_field_tecplot(const Qpoint_func& output_variables, std::string name, int n_sample = 20,
                               bool edges = false, bool qpoints = false, bool interior = true);
  /*! if a `Qpoint_func` is not specified, all the state variables
   * and the artificial viscosity coefficient will be output
   */
  void visualize_field_tecplot(std::string name, int n_sample = 20,
                               bool edges = false, bool qpoints = false, bool interior = true);
  /*! write a visualization file describing all surfaces where a particular boundary condition has been enforced.
   * only does state variables because usually that's what you want and I'm lazy
   */
  void visualize_surface_tecplot(int bc_sn, std::string name, int n_sample = 20);
  #endif
  #if HEXED_USE_OTTER
  void visualize_edges_otter(otter::plot&, Eigen::Matrix<double, 1, Eigen::Dynamic> color = otter::colors::css4["white"], int n_sample = 21);
  /*!
   * plot the surface with optional color mapping. plotting takes whatever form is appropriate for dimensionality
   * note: color_by must be a scalar
   * if either element of `bounds` is NaN, will substitute min & max of variable in domain
   */
  void visualize_surface_otter(otter::plot&, int bc_sn, const otter::colormap& = otter::const_colormap(otter::colors::css4["darkgrey"]),
                               const Qpoint_func& color_by = Pressure(), std::array<double, 2> bounds = {std::nan(""), std::nan("")}, bool transparent = false, int n_sample = 21, double tol = 1e-3);
  /*!
   * plot the flow field in the most appropriate way for the dimensionality
   * includes contour lines/surfaces, and for 2d also colors flowfield
   * note that this can also be used to plot slices if you contour by a `Linear` object
   */
  void visualize_field_otter(otter::plot&,
                             const Qpoint_func& contour = Pressure(),
                             int n_contour = 3,
                             std::array<double, 2> contour_bounds = {std::nan(""), std::nan("")},
                             const Qpoint_func& color_by = Pressure(),
                             std::array<double, 2> color_bounds = {std::nan(""), std::nan("")},
                             const otter::colormap& cmap_contour = otter::const_colormap(Eigen::Vector4d{1., 1., 1., .1}),
                             const otter::colormap& cmap_field   = otter::plasma,
                             bool transparent = true, bool show_color = true,
                             int n_div = 10, double tol = 1e-3);
  /*!
   * convenience function to plot a slice
   * wrapper for `visualize_field_otter`
   */
  inline void visualize_slice_otter(otter::plot& plt, int i_dim, double location,
                                    const Qpoint_func& color_by = Pressure(),
                                    std::array<double, 2> bounds = {std::nan(""), std::nan("")},
                                    const otter::colormap& cmap = otter::plasma,
                                    int n_div = 10)
  {
    visualize_field_otter(plt, hexed::Linear(Eigen::VectorXd::Unit(params.n_dim, i_dim)),
                          1, {location, location},
                          color_by, bounds, cmap, cmap, false, n_div);
  }
  #endif
  //!\}
};

}
#endif
