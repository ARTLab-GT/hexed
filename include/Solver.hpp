#ifndef HEXED_SOLVER_HPP_
#define HEXED_SOLVER_HPP_

#include "Solver_interface.hpp"
#include "Accessible_mesh.hpp"
#include "kernel_factory.hpp"
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
class Solver : public Solver_interface
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

  Mesh& mesh() override;
  Storage_params storage_params() override;
  void relax_vertices(double factor = .5) override;
  void snap_vertices() override;
  void snap_faces() override;
  void calc_jacobian() override;
  void initialize(const Spacetime_func&) override;
  void set_art_visc_off() override;
  void set_art_visc_constant(double) override;
  void set_art_visc_row_size(int) override;
  void set_fix_admissibility(bool) override;
  void set_resolution_badness(const Element_func& func) override;
  void set_res_bad_surface_rep(int bc_sn) override;
  void synch_extruded_res_bad() override;
  void update(double dt = 0.8, bool cfl_driven = true) override;
  bool is_admissible() override;
  void set_art_visc_smoothness(double advect_length) override;
  Iteration_status iteration_status() override;
  void reset_counters() override;
  std::vector<double> sample(int ref_level, bool is_deformed, int serial_n, int i_qpoint, const Qpoint_func&) override;
  std::vector<double> sample(int ref_level, bool is_deformed, int serial_n, const Element_func&) override;
  const Stopwatch_tree& stopwatch_tree() override;
  std::vector<double> integral_field(const Qpoint_func& integrand) override;
  std::vector<double> integral_surface(const Boundary_func& integrand, int bc_sn) override;
  std::vector<std::array<double, 2>> bounds_field(const Qpoint_func&, int n_sample = 20) override;
  #if HEXED_USE_TECPLOT
  void visualize_field_tecplot(const Qpoint_func& output_variables, std::string name, int n_sample = 20,
                               bool edges = false, bool qpoints = false, bool interior = true) override;
  void visualize_field_tecplot(std::string name, int n_sample = 20,
                               bool edges = false, bool qpoints = false, bool interior = true) override;
  void visualize_surface_tecplot(int bc_sn, const Boundary_func&, std::string name, int n_sample = 20) override;
  void visualize_surface_tecplot(int bc_sn, std::string name, int n_sample = 20) override;
  void vis_cart_surf_tecplot(int bc_sn, std::string name, const Boundary_func& func = Resolution_badness()) override;
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
};

}
#endif
