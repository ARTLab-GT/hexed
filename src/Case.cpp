#include <filesystem>
#include <Case.hpp>

namespace hexed
{

Solver& Case::_solver()
{
  HEXED_ASSERT(_solver_ptr, "`Solver` object does not exist");
  return *_solver_ptr;
}

std::optional<        int> Case::_vari(std::string name) {return _inter.variables->lookup<        int>(name);}
std::optional<     double> Case::_vard(std::string name) {return _inter.variables->lookup<     double>(name);}
std::optional<std::string> Case::_vars(std::string name) {return _inter.variables->lookup<std::string>(name);}

Mat<> Case::_get_vector(std::string name, int size)
{
  Mat<> vec(size);
  for (int i = 0; i < size; ++i) {
    HEXED_ASSERT(_vard(name + std::to_string(i)), "must specify all components of `" + name + "` or none");
    vec(i) = *_vard(name + std::to_string(i));
  }
  return vec;
}
void Case::_set_vector(std::string name, Mat<> vec)
{
  for (int i = 0; i < int(vec.size()); ++i) {
    _inter.variables->assign<double>(name + std::to_string(i), vec(i));
  }
}

Case::Case(std::string input_file)
{
  // create custom Heisenberg variables
  _inter.variables->create<int>("create_solver", new Namespace::Heisenberg<int>([this]() {
    auto n_dim = _inter.variables->lookup<int>("n_dim");
    HEXED_ASSERT(n_dim && n_dim.value() > 0 && n_dim.value() <= 3,
                 "`n_dim` must be defined as an integer in [1, 3]");
    auto row_size = _inter.variables->lookup<int>("row_size");
    HEXED_ASSERT(row_size.value() >= 2 && row_size.value() <= config::max_row_size,
                 format_str(300, "`row_size` must be between 2 and %i", config::max_row_size));
    // compute freestream
    Mat<> freestream(*n_dim + 2);
    double heat_rat = 1.4;
    if (_vard("freestream0")) freestream = _get_vector("freestream", *n_dim + 2);
    else {
      HEXED_ASSERT(_vard("density").has_value() + _vard("pressure").has_value() + _vard("temperature").has_value() == 2,
                   "exactly two of density, pressure, and temperature must be specified");
      if (_vard("density")) {
        freestream(*n_dim) = *_vard("density");
        if (_vard("pressure")) _inter.variables->assign<double>("temperature", *_vard("pressure")/(constants::specific_gas_air**_vard("density")));
        else _inter.variables->assign<double>("pressure", *_vard("density")*constants::specific_gas_air**_vard("temperature"));
      } else _inter.variables->assign<double>("density", *_vard("pressure")/constants::specific_gas_air**_vard("temperature"));
      HEXED_ASSERT(_vard("velocity0").has_value() + _vard("speed").has_value() + _vard("mach").has_value() == 1,
                   "exactly one of velosity, speed, and Mach number must be specified");
      Mat<> veloc;
      if (_vard("velocity0")) veloc = _get_vector("velocity", *n_dim);
      else {
        Mat<> direction;
        if (_vard("direction0")) direction = _get_vector("direction", *n_dim).normalized();
        else {
          direction.setUnit(*n_dim, 0);
          if (*n_dim == 2) {
            direction = Eigen::Rotation2D<double>(_vard("attack").value())*direction;
          }
          if (*n_dim == 3) {
            direction = Eigen::AngleAxis<double>(-_vard("attack"  ).value(), Eigen::Vector3d::Unit(1))*direction;
            direction = Eigen::AngleAxis<double>( _vard("sideslip").value(), Eigen::Vector3d::Unit(2))*direction;
          }
          _set_vector("direction", direction);
          _inter.variables->assign<double>("sound", std::sqrt(heat_rat*constants::specific_gas_air**_vard("temperature")));
          if (_vard("speed")) _inter.variables->assign<double>("mach", *_vard("speed")/ *_vard("sound"));
          else _inter.variables->assign<double>("speed", *_vard("mach")**_vard("sound"));
        }
        veloc = *_vard("speed")*direction;
        _set_vector("velocity", veloc);
      }
      freestream(Eigen::seqN(0, *n_dim)) = *_vard("density")*veloc;
      freestream(*n_dim) = *_vard("density");
      freestream(*n_dim + 1) = *_vard("pressure")/(heat_rat - 1) + .5**_vard("density")*veloc.squaredNorm();
      _set_vector("freestream", freestream);
    }
    // create solver
    Mat<dyn, dyn> mesh_corners(*n_dim, 2);
    std::vector<Flow_bc*> bcs;
    for (int i_dim = 0; i_dim < *n_dim; ++i_dim) {
      for (int sign = 0; sign < 2; ++sign) {
        std::string index = format_str(50, "%i%i", i_dim, sign);
        mesh_corners(i_dim, sign) = _inter.variables->lookup<double>("mesh_corner" + index).value();
        std::string bc_name = _inter.variables->lookup<std::string>("extremal_bc" + index).value();
        if (bc_name == "characteristic") bcs.push_back(new Riemann_invariants(freestream));
        else if (bc_name == "nonpenetration") bcs.push_back(new Nonpenetration);
        else HEXED_ASSERT(false, format_str(1000, "unrecognized boundary condition type `%s`", bc_name));
      }
    }
    HEXED_ASSERT((mesh_corners(all, 1) - mesh_corners(all, 0)).minCoeff() > 0, "all mesh dimensions must be positive!");
    double root_sz = (mesh_corners(all, 1) - mesh_corners(all, 0)).maxCoeff();
    _solver_ptr.reset(new Solver(*n_dim, *row_size, root_sz));
    _solver().mesh().add_tree(bcs, mesh_corners(all, 0));
    return 0;
  }));

  _inter.variables->create<int>("init_refinement", new Namespace::Heisenberg<int>([this]() {
    for (int i = 0; i < _inter.variables->lookup<int>("init_ref_level"); ++i) _solver().mesh().update();
    return 0;
  }));

  _inter.variables->create<int>("init_state", new Namespace::Heisenberg<int>([this]() {
    std::string init_cond = _vars("init_condition").value();
    auto freestream = _get_vector("freestream", _vari("n_dim").value() + 2);
    if (init_cond == "freestream") {
      _solver().initialize(Constant_func({freestream.begin(), freestream.end()}));
    } else if (init_cond == "vortex") {
      Isentropic_vortex vortex({freestream.begin(), freestream.end()});
      vortex.max_nondim_veloc = 0.3;
      vortex.argmax_radius = 0.1;
      _solver().initialize(vortex);
    } else HEXED_ASSERT(false, format_str(1000, "unrecognized initial condition type `%s`", init_cond.c_str()));
    return 0;
  }));

  _inter.variables->create<int>("visualize", new Namespace::Heisenberg<int>([this]() {
    std::filesystem::path file_name(_inter.variables->lookup<std::string>("vis_file_name").value());
    if (!std::filesystem::exists(file_name.parent_path())) {
      std::filesystem::create_directory(file_name.parent_path());
    }
    _solver().visualize_field_tecplot(file_name.string());
    return 0;
  }));

  // load HIL code for the Case _interface
  _inter.exec(format_str(1000, "$read \"%s/include/Case.hil\"", config::root_dir));
  // execute input file
  _inter.exec(format_str(1000, "$read \"%s\"", input_file.c_str()));
}

}
