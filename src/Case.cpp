#include <filesystem>
#include <Case.hpp>
#include <Simplex_geom.hpp>

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
    Mat<dyn, dyn> mesh_bounds(*n_dim, 2);
    std::vector<Flow_bc*> bcs;
    for (int i_dim = 0; i_dim < *n_dim; ++i_dim) {
      for (int sign = 0; sign < 2; ++sign) {
        std::string index = format_str(50, "%i%i", i_dim, sign);
        mesh_bounds(i_dim, sign) = _inter.variables->lookup<double>("mesh_bound" + index).value();
        std::string bc_name = _inter.variables->lookup<std::string>("extremal_bc" + index).value();
        if (bc_name == "characteristic") bcs.push_back(new Freestream(freestream));
        else if (bc_name == "nonpenetration") bcs.push_back(new Nonpenetration);
        else HEXED_ASSERT(false, format_str(1000, "unrecognized boundary condition type `%s`", bc_name));
      }
    }
    HEXED_ASSERT((mesh_bounds(all, 1) - mesh_bounds(all, 0)).minCoeff() > 0, "all mesh dimensions must be positive!");
    double root_sz = (mesh_bounds(all, 1) - mesh_bounds(all, 0)).maxCoeff();
    _solver_ptr.reset(new Solver(*n_dim, *row_size, root_sz));
    _solver().mesh().add_tree(bcs, mesh_bounds(all, 0));
    return 0;
  }));

  _inter.variables->create<int>("init_refinement", new Namespace::Heisenberg<int>([this]() {
    for (int i = 0; i < _inter.variables->lookup<int>("init_ref_level"); ++i) _solver().mesh().update();
    _solver().calc_jacobian();
    return 0;
  }));

  _inter.variables->create<int>("add_geom", new Namespace::Heisenberg<int>([this]() {
    int nd = _vari("n_dim").value();
    std::vector<Surface_geom*> geoms;
    for (int i_geom = 0;; ++i_geom) {
      std::string prefix = "geom_" + std::to_string(i_geom);
      auto geom_type = _vars(prefix + "_type");
      if (!geom_type) break;
      if (geom_type.value() == "file") {
        HEXED_ASSERT(false, "geometry from file names is not yet supported");
      } else if (geom_type.value() == "parametric") {
        HEXED_ASSERT(nd == 2,
          "Parametric geometry is currently only supported for `n_dim == 2`. "
          "If you would like this for other dimensionality, lmk.");
        auto def = _vars(prefix + "_definition");
        auto n = _vari(prefix + "_n");
        HEXED_ASSERT(def && n, format_str(200, "specifying geometry %i requires `geom_%i_n` and `geom_%i_definition` to be specified", i_geom, i_geom, i_geom));
        Mat<dyn, dyn> points (nd, *n + 1);
        for (int i_point = 0; i_point <= *n; ++i_point) {
          auto sub = _inter.make_sub();
          sub.exec(format_str(1000, "param = (0. + %i)/%i; $equations", i_point, *n));
          for (int i_dim = 0; i_dim < nd; ++i_dim) {
            auto val = sub.variables->lookup<double>("x_" + std::to_string(i_dim));
            HEXED_ASSERT(val, "parametric geometry equations failed to set all coordinates");
            points(i_dim, i_point) = val.value();
          }
        }
        geoms.push_back(new Simplex_geom<2>(segments(points)));
      } else HEXED_ASSERT(false, format_str(1000, "geometry type `%s` not recognized", geom_type.value()));
    }
    if (!geoms.empty()) _solver().mesh().set_surface(new Compound_geom(geoms), new Nonpenetration, _get_vector("flood_fill_start", nd));
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

  _inter.variables->create<int>("update", new Namespace::Heisenberg<int>([this]() {
    _solver().update(_vard("max_safety").value(), _vard("max_time_step").value());
    return 0;
  }));
  _inter.variables->create<std::string>("header", new Namespace::Heisenberg<std::string>([this]() {
    return _solver().iteration_status().header();
  }));
  _inter.variables->create<std::string>("report", new Namespace::Heisenberg<std::string>([this]() {
    return _solver().iteration_status().report();
  }));

  // load HIL code for the Case _interface
  _inter.exec(format_str(1000, "$read {%s/include/Case.hil}", config::root_dir));
  // execute input file
  _inter.exec(format_str(1000, "$read {%s}", input_file.c_str()));
}

}
