#ifndef HEXED_SOLVER_INTERFACE_HPP_
#define HEXED_SOLVER_INTERFACE_HPP_

#include "config.hpp"
#include "Mesh_interface.hpp"
#include "Gauss_legendre.hpp"
#include "Spacetime_func.hpp"
#include "Element_func.hpp"
#include "Iteration_status.hpp"
#include "Stopwatch_tree.hpp"
#include "Transport_model.hpp"

namespace hexed
{

class Solver_interface
{
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

  virtual ~Solver_interface() = default;

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
  virtual Mesh_interface& mesh() = 0;
  virtual Storage_params storage_params() = 0;
  //! moves all vertices to the mean of the current position and the mean of the neighbors' positions
  virtual void relax_vertices(double factor = .5) = 0; //!< \deprecated use `Mesh::relax`
  /*! \brief apply `Mesh_bc`s
   * This is where hanging vertices are snapped to their coarse faces.
   * \note if some elements participate in multiple BCs, then snapping may not satisfy all exactly.
   *       However, if performed multiple times, it should converge = 0;
   */
  virtual void snap_vertices() = 0; //!< \deprecated tree meshing does this automatically
  //! warps the boundary elements such that the element faces coincide with the boundary at their quadrature points.
  virtual void snap_faces() = 0;
  /*!
   * compute the Jacobian of all elements based on the current position of the vertices
   * and value of any face warping.
   * Mesh topology must be valid (no duplicate or missing connections) before calling this function.
   */
  virtual void calc_jacobian() = 0;
  //! set the flow state
  virtual void initialize(const Spacetime_func&) = 0;
  virtual void set_art_visc_off() = 0; //!< turns off artificial viscosity
  virtual void set_art_visc_constant(double) = 0; //!< turns on artificial viscosity and initializes coefficient to a uniform value
  virtual void set_art_visc_row_size(int) = 0; //!< modify the polynomial order of smoothness-based artificial viscosity (must be <= row size of discretization (which is the default))
  virtual void set_fix_admissibility(bool) = 0; //!< turns on/off the thermodynamic admissibility-preserving scheme (increases robustness at some computational overhead)
  /*! \brief set `Element::resolution_badness` for each element according to `func`.
   * \details Resoltuion badness can be evaluated via `sample(ref_level, is_deformed, serial_n, Resolution_badness())`.
   * This function does some additional work to enforce some conditions on the resolution badness of neighboring elements.
   * Thus, use this function rather than just `sample(ref_level, is_deformed, serial_n, func)` directly.
   */
  virtual void set_resolution_badness(const Element_func& func) = 0;
  /*! \brief Set resolution badness based on surface representation quality.
   * \details For all deformed elements contacting the boundary specified by `bc_sn`, computes the
   * unit surface normals at the faces and compares to neighboring elements.
   * `Element::resolution_badness` is set to the total difference between unit normals with all neighbors,
   * where in 3D the total on each edge is computed by Gaussian quadrature in reference space.
   * \note Connections between deformed elements and Cartesian elements are ignored.
   * Fully Cartesian connections trivially contribute zero to this metric of resolution badness.
   */
  virtual void set_res_bad_surface_rep(int bc_sn) = 0;
  //! set resolution badness of each element to be at least the maximum badness of any elements extruded from it
  virtual void synch_extruded_res_bad() = 0;
  //!\}

  //! \name time marching
  //!\{
  /*!
   * March the simulation forward by a time step equal to `time_step` or
   * `safety_factor` times the estimated maximum stable time step, whichever is smaller.
   * \note For local time stepping, `time_step` is equivalent to the CFL number.
   * Also, `safety_factor` is __not__ the same as the CFL number
   * (it is scaled by the max allowable CFL for the chosen DG scheme which is often O(1e-2)).
   */
  virtual void update(double safety_factor = 0.7, double time_step = std::numeric_limits<double>::max()) = 0;
  //! runs `update` `n` times, in case you're calling it through a slow interface _cough_ python _cough_
  virtual void update_n(int n, double safety_factor = 0.7, double time_step = std::numeric_limits<double>::max()) = 0;
  virtual bool is_admissible() = 0; //!< check whether flowfield is admissible (e.g. density and energy are positive)
  virtual void set_art_visc_smoothness(double advect_length) = 0; //!< updates the aritificial viscosity coefficient based on smoothness of the flow variables
  /*! \brief an object providing all available information about the status of the time marching iteration.
   * \details The `Iteration_status::start_time` member will refer to when the `Solver` object was created
   * (specifically at the start of the `Solver::Solver` body).
   */
  virtual Iteration_status iteration_status() = 0;
  /*!
   * reset any variables in `iteration_status()` that count something since the last call to `reset_counters()`
   * (e.g. number of artificial viscosity iterations)
   */
  virtual void reset_counters() = 0;
  //!\}

  //! \name output
  //! functions that compute some form of output data
  //!\{
  virtual std::vector<double> sample(int ref_level, bool is_deformed, int serial_n, int i_qpoint, const Qpoint_func&) = 0; //!< evaluate arbitrary functions at arbitrary locations
  virtual std::vector<double> sample(int ref_level, bool is_deformed, int serial_n, const Element_func&) = 0; //!< \overload
  //! obtain performance data
  virtual const Stopwatch_tree& stopwatch_tree() = 0;
  //! compute an integral over the entire flow field at the current time
  virtual std::vector<double> integral_field(const Qpoint_func& integrand) = 0;
  //! compute an integral over all surfaces where a particular boundary condition has been enforced
  virtual std::vector<double> integral_surface(const Boundary_func& integrand, int bc_sn) = 0;
  /*! compute the min and max of variables over entire flow field. layout: `{{var0_min, var0_max}, {var1_min, var1_max}, ...}`
   * bounds are approximated by uniformly sampling a block `n_sample`-on-a-side in each element
   */
  virtual std::vector<std::array<double, 2>> bounds_field(const Qpoint_func&, int n_sample = 20) = 0;
  #if HEXED_USE_TECPLOT
  //! write a visualization file describing the entire flow field (but not identifying surfaces)
  virtual void visualize_field_tecplot(const Qpoint_func& output_variables, std::string name, int n_sample = 20,
                               bool edges = false, bool qpoints = false, bool interior = true) = 0;
  /*! if a `Qpoint_func` is not specified, all the state variables
   * and the artificial viscosity coefficient will be output
   */
  virtual void visualize_field_tecplot(std::string name, int n_sample = 20,
                               bool edges = false, bool qpoints = false, bool interior = true) = 0;
  /*! write a visualization file describing all surfaces where a particular boundary condition has been enforced.
   */
  virtual void visualize_surface_tecplot(int bc_sn, const Boundary_func&, std::string name, int n_sample = 20) = 0;
  /*! if a `Boundary_func` is not specified, all the state variables,
   * the surface normal, the shear stress, and the heat flux will be output.
   */
  virtual void visualize_surface_tecplot(int bc_sn, std::string name, int n_sample = 20) = 0;
  //! visualize the Cartesian surface which theoretically exists after element deletion but before any vertex snapping
  virtual void vis_cart_surf_tecplot(int bc_sn, std::string name, const Boundary_func& func = Resolution_badness()) = 0;
  #endif
};

std::unique_ptr<Solver_interface> make_solver(int n_dim, int row_size, double root_mesh_size, bool local_time_stepping = false,
                                              Transport_model viscosity_model = inviscid, Transport_model thermal_conductivity_model = inviscid);

}
#endif
