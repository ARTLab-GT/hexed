#ifndef HEXED_SOLVER_HPP_
#define HEXED_SOLVER_HPP_

#include "config.hpp"
#include "Gauss_legendre.hpp"
#include "Spacetime_func.hpp"
#include "Element_func.hpp"
#include "Iteration_status.hpp"
#include "Stopwatch_tree.hpp"
#include "Transport_model.hpp"
#include "Mesh.hpp"
#include "Accessible_mesh.hpp"
#include "Namespace.hpp"
#include "Printer.hpp"
#include "Linear_equation.hpp"
#include "Visualizer.hpp"
#include "kernels.hpp"

namespace hexed
{

/*! \brief The main class that basically runs everything.
 * \details If you want to run a simulation with hexed through the C++ API, your workflow should be roughly the following:
 * -# construct a `Solver` object
 * -# interact with the `mesh()` object to build the mesh topology and/or snap vertices
 * -# call `calc_jacobian()` to initialize internal parameters based on the mesh
 * -# call `initialize()` to initialize the flow state
 * -# call `update()` repeatedly to progress the simulation
 * -# call the some of the functions in output section to get the data you want from the simulation
 */
class Solver
{
  Storage_params params;
  std::unique_ptr<Accessible_mesh> acc_mesh;
  Gauss_legendre basis;
  Iteration_status status;
  Stopwatch_tree stopwatch;
  bool use_art_visc;
  bool fix_admis;
  int av_rs;
  Transport_model visc;
  Transport_model therm_cond;
  int last_fix_vis_iter = std::numeric_limits<int>::min();
  std::shared_ptr<Namespace> _namespace;
  std::shared_ptr<Printer> _printer;
  bool _implicit;

  Kernel_mesh _kernel_mesh();
  void share_vertex_data(std::function<double&(Element&, int i_vertex)>, std::function<double(Mat<>)>);
  void share_vertex_data(std::function<double(Element&, int i_vertex)> get, std::function<double&(Element&, int i_vertex)> set, std::function<double(Mat<>)>);
  bool fix_admissibility(double stability_ratio);
  void apply_state_bcs();
  void apply_flux_bcs();
  void apply_avc_diff_bcs();
  void apply_avc_diff_flux_bcs();
  void apply_fta_flux_bcs();
  void diffuse_art_visc(double diff_time);
  void fta(double dt, int i_stage);
  bool use_ldg();
  double max_dt(double max_safety_conv, double max_safety_diff);
  std::unique_ptr<Visualizer> _visualizer(std::string format, std::string name, const Output_data& output_variables, int n_dim_topo);

  //! \brief linearizes the steady state equations by finite difference
  class Linearized : public Linear_equation
  {
    Solver& _solver;
    Mat<> _ref_state;
    Mat<> _weights;
    public:
    //! \brief the first `storage_start` vectors are reserved for internal use
    //! \details (0, 1, and 2 are the working data for the Spatial kernel and 3 is the linearization point)
    static const int storage_start = 4;
    static const int n_storage = 30; //!< \brief number of vectors available for algorithms to work with
    double finite_diff = 1e-3; //!< \brief \f$ x \f$ vectors will be scaled by this amount to improve linearity
    Linearized(Solver&); //!< \brief sets linearization point to current working state
    int n_vecs() override;
    void scale(int output, int input, double scalar) override;
    void add(int output, double coef0, int vec0, double coef1, int vec1) override;
    double inner(int input0, int input1) override;
    /*! \brief applies linearized operator
     * \details defined as \f$ A x = \frac{1}{\epsilon}(\nabla F(u) - \nabla F(u + \epsilon x)) \f$
     * where \f$ u \f$ is the linearization point (the state to linearize about), \f$ \epsilon \f$ is given by `Linearized::finite_diff`,
     * and \f$ \nabla F \f$ is of course the residual of the nonlinear steady-state equations.
     */
    void matvec(int output, int input) override;
  };

  public:
  /*!
   * \param n_dim number of dimensions
   * \param row_size row size of the basis (see \ref Terminology)
   * \param root_mesh_size sets the value of `Mesh::root_mesh_size()`
   * \param local_time_stepping whether to use local or global time stepping
   * \param viscosity_model determines whether the flow has viscosity (natural, not artificial) and if so, how it depends on temperature
   * \param thermal_conductivity_model determines whether the flow has thermal conductivity and if so, how it depends on temperature
   * \param space `Namespace` containing any user-defined parameters affecting the behavior of the solver.
   *        If no namespace is provided, a new blank namespace is creqated.
   *        Any optional parameters which are not found in the namespace shall be created with their default values.
   * \param printer what to do with any information the solver wants to print for the user to see
   * \param implicit if `true`, allocate storage for solving with an implicit method (experimental feature -- not ready for production use)
   * \details If `viscosity_model` and `thermal_conductivity_model` are both `inviscid` _and_ you don't turn on artificial viscosity,
   * you will be solving the pure inviscid flow equations.
   * Otherwise, you will be solving the viscous flow equations using the LDG scheme,
   * potentially with some of the diffusion coefficients
   * (artificial viscosity, natural viscosity, thermal conductivity) set to zero.
   */
  Solver(int n_dim, int row_size, double root_mesh_size, bool local_time_stepping = false,
         Transport_model viscosity_model = inviscid, Transport_model thermal_conductivity_model = inviscid,
         std::shared_ptr<Namespace> space = std::make_shared<Namespace>(), std::shared_ptr<Printer> printer = std::make_shared<Stream_printer>(),
         bool implicit = false);

  //! \name setup
  //!\{
  Namespace& nspace(); //! \brief Reference to the namespace which can be used to edit user-defined parameters
  /*! \brief fetch the `Mesh`.
   * \details An object the user can use to build the mesh.
   * Note that whenever elements are added, the flow state, and Jacobian are uninitialized,
   * the time step scale is uniformly 1,
   * and the mesh quality may be poor.
   * The functions below must be used to complete the setup
   * before any flow calculation can begin.
   */
  Mesh& mesh();
  /*! \brief reads mesh from file
   * \details Wipes old mesh and flow state.
   * File must be in the native (HDF5-based) mesh format, which you can generate from a previous simulation with `Mesh::write`.
   * New mesh must match the `Storage_params` of the current one, but the root mesh size will be replaced with that of the new mesh.
   * You still have to `initialize` (even if you already did before reading the new mesh),
   * but you don't have to `calc_jacobian` unless you further modify the mesh.
   */
  void read_mesh(std::string file_name, std::vector<Flow_bc*> extremal_bcs, Surface_geom* = nullptr, Flow_bc* surface_bc = nullptr);
  Storage_params storage_params();
  //! warps the boundary elements such that the element faces coincide with the boundary at their quadrature points.
  void snap_faces();
  /*! \brief compute the Jacobian of all elements based on the current position of the vertices and value of any face warping.
   * \details Mesh topology must be valid (no duplicate or missing connections) before calling this function.
   * \param `snap_faces` if `true`, this function will go ahead and perform face snapping for you
   */
  void calc_jacobian(bool snap_faces = true);
  //! set the flow state
  void initialize(const Spacetime_func&);
  void set_art_visc_off(); //!< turns off artificial viscosity
  void set_art_visc_constant(double); //!< turns on artificial viscosity and initializes coefficient to a uniform value
  void set_art_visc_row_size(int); //!< modify the polynomial order of smoothness-based artificial viscosity (must be <= row size of discretization (which is the default))
  void set_fix_admissibility(bool); //!< turns on/off the thermodynamic admissibility-preserving scheme (increases robustness at some computational overhead)
  /*! \brief set `Element::uncertainty` for each element according to `func`.
   * \details Uncertainty metric can be evaluated via `sample(ref_level, is_deformed, serial_n, Uncertainty())`.
   * This function does some additional work to enforce some conditions on the uncertainty of neighboring elements.
   * Thus, use this function rather than just `sample(ref_level, is_deformed, serial_n, func)` directly.
   */
  void set_uncertainty(const Element_func& func);
  /*! \brief Set uncertainty metric based on surface representation quality.
   * \details For all deformed elements contacting the boundary specified by `bc_sn`, computes the
   * unit surface normals at the faces and compares to neighboring elements.
   * `Element::uncertainty` is set to the total difference between unit normals with all neighbors,
   * where in 3D the total on each edge is computed by Gaussian quadrature in reference space.
   * \note Connections between deformed elements and Cartesian elements are ignored.
   * Fully Cartesian connections trivially contribute zero to this metric of uncertainty.
   */
  void set_uncert_surface_rep(int bc_sn);
  //! set uncertainty of each element to be at least the maximum uncertainty of any elements extruded from it
  void synch_extruded_uncert();
  //!\}

  //! \name time marching
  //!\{
  /*!
   * March the simulation forward by a time step equal to `time_step` or
   * `max_safety` times the estimated maximum stable time step, whichever is smaller.
   * Also, the safety factor is __not__ the same as the CFL number
   * (it is scaled by the max allowable CFL for the chosen DG scheme which is often O(1e-2)).
   */
  void update();
  void update_implicit(); //!< \brief (experimental) performs an implicit time step \warning Experimental! Interesting for reasearch, not effective in practice (yet, anyway).
  void compute_residual();
  bool is_admissible(); //!< \brief check whether flowfield is admissible (e.g. density and energy are positive)
  void update_art_visc_smoothness(double advect_length); //!< \brief updates the aritificial viscosity coefficient based on smoothness of the flow variables
  /*! \brief (experimental) sets artificial viscosity based on elementwise smoothness
   * \details Based on the work of Persson et al. on artificial viscosity with elementwise smoothness indicators.
   * Implemented primarily for evaluating the difference between grid-independent artificial viscosity and conventional methods.
   * \warning Not recommended for use in practice.
   * Use `Solver::update_art_visc_smoothness`, which was developed to supplant this type of approach.
   */
  void update_art_visc_elwise(double width, bool pde_based = false);
  /*! \brief an object providing all available information about the status of the time marching iteration.
   * \details The `Iteration_status::start_time` member will refer to when the `Solver` object was created
   * (specifically at the start of the `Solver::Solver` body).
   */
  void set_art_visc_admis(); //!< \brief set the Laplacian artificial viscosity to a minimal value that will encourage thermodynamic admissibility
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
  std::vector<double> sample(int ref_level, bool is_deformed, int serial_n, const Element_func&); //!< \overload
  //! obtain performance data
  const Stopwatch_tree& stopwatch_tree();
  //! compute an integral over the entire flow field at the current time
  std::vector<double> integral_field(const Qpoint_func& integrand);
  //! compute an integral over all surfaces where a particular boundary condition has been enforced
  std::vector<double> integral_surface(const Boundary_func& integrand, int bc_sn);
  /*! compute the min and max of variables over entire flow field. layout: `{{var0_min, var0_max}, {var1_min, var1_max}, ...}`
   * bounds are approximated by uniformly sampling a block `n_sample`-on-a-side in each element
   */
  std::vector<std::array<double, 2>> bounds_field(const Qpoint_func&, int n_sample = 20);

  /*! \brief write a visualization file describing the entire flow field (but not identifying surfaces)
   * \param format Which format to write the visualization file in. Accepted values are `"xdmf"` and `"tecplot"`
   * \param name name of file to write (not including extension)
   * \param output_variables what variables to write
   * \param n_sample each element will contain an `n_sample` by `n_sample` array of uniformly-spaced sample points
   * \param wireframe if `true`, visualize the mesh edges as a wireframe instead of the filled surface/solid
   */
  void visualize_field(std::string format, std::string name, const Qpoint_func& output_variables, int n_sample = 10, bool wireframe = false);
  //! \brief write a visualization file describing all surfaces where a particular boundary condition has been enforced.
  void visualize_surface(std::string format, std::string name, int bc_sn, const Boundary_func&, int n_sample = 10, bool wireframe = false);
  //! \brief visualize the Cartesian surface which theoretically exists after element deletion but before any vertex snapping
  void vis_cart_surf(std::string format, std::string name, int bc_sn, const Boundary_func& func = Uncertainty());
  //! \brief visualize the local time step constraints imposed by convection and diffusion, respectively
  //! \warning This function overwrites the reference state, which will invalidate any residual evaluation until `update` is called again.
  void vis_lts_constraints(std::string format, std::string name, int n_sample = 10);
  //!\}
};

}
#endif
