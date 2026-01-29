#include <filesystem>
#include <cctype>
#include <hexed/Case.hpp>
#include <hexed/Simplex_geom.hpp>
#include <hexed/read_csv.hpp>
#include <hexed/standard_atmosphere.hpp>
#include <hexed/vis_variables.hpp>
#include <hexed/Csv.hpp>
#include <hexed/brep.hpp>
#include <hexed/Printer.hpp>
#include <hexed/Tree_curve_geom.hpp>

namespace hexed {

const double heat_rat = 1.4;

Solver& Case::_solver() {
  HEXED_ASSERT(_solver_ptr, "`Solver` object does not exist", assert::User_error);
  return *_solver_ptr;
}

std::string strip_trailing_digits(std::string s) {
  while (std::isdigit(s.back())) s.pop_back();
  return s;
}

        int Case::_vari(std::string name) {return _inter.variables->get<        int>(name);}
     double Case::_vard(std::string name) {return _inter.variables->get<     double>(name);}
std::string Case::_vars(std::string name) {return _inter.variables->get<std::string>(name);}

Mat<> Case::_get_vector(std::string name, int size) {
  Mat<> vec(size);
  for (int i = 0; i < size; ++i) {
    HEXED_ASSERT(_inter.variables->lookup<double>(name + std::to_string(i)),
                 "must specify all components of `" + name + "` or none", assert::User_error)
    vec(i) = _vard(name + std::to_string(i));
  }
  return vec;
}
void Case::_set_vector(std::string name, Mat<> vec) {
  for (int i = 0; i < int(vec.size()); ++i) {
    _inter.variables->assign<double>(name + std::to_string(i), vec(i));
  }
}

std::shared_ptr<Flow_bc> Case::_make_bc(std::string name) {
  std::shared_ptr<Flow_bc> bc;
  Mat<> freestream = _get_vector("freestream", _vari("n_var"));
  if      (name == "characteristic") bc = std::make_shared<Riemann_invariants>(freestream);
  else if (name == "freestream") bc = std::make_shared<Freestream>(freestream);
  else if (name == "pressure_outflow") bc = std::make_shared<Pressure_outflow>(_vard("freestream_pressure"));
  else if (name == "outflow") bc = std::make_shared<Outflow>();
  else if (name == "nonpenetration") bc = std::make_shared<Nonpenetration>();
  else if (name == "no_slip") {
    auto sub = _inter.make_sub();
    sub.exec("$thermal_bc");
    std::shared_ptr<Thermal_bc> thermal;
    if (sub.variables->exists("heat_flux")) { // note not recursive
      thermal = std::make_shared<Prescribed_heat_flux>(sub.variables->get<double>("heat_flux"));
    } else if (sub.variables->exists("emissivity") || sub.variables->exists("heat_transfer_coef")) {
      auto equilibrium = std::make_shared<Thermal_equilibrium>();
      HEXED_ASSERT(sub.variables->exists("heat_transfer_coef") == sub.variables->exists("temperature"),
        "must specify both surface heat_transfer_coef and temperature or neither", assert::User_error);
      if (sub.variables->exists("emissivity")) equilibrium->emissivity = sub.variables->get<double>("emissivity");
      if (sub.variables->exists("heat_transfer_coef")) {
        equilibrium->heat_transfer_coef = sub.variables->get<double>("heat_transfer_coef");
        equilibrium->temperature = sub.variables->get<double>("temperature");
        equilibrium->heat_rat = heat_rat;
      }
      thermal = equilibrium;
    } else if (sub.variables->exists("internal_energy")) { // note not recursive
      thermal = std::make_shared<Prescribed_energy>(sub.variables->get<double>("internal_energy"));
    } else if (sub.variables->exists("temperature")) { // note not recursive
      double energy = sub.variables->get<double>("temperature")*constants::specific_gas_air/(heat_rat - 1.);
      thermal = std::make_shared<Prescribed_energy>(energy);
    }
    HEXED_ASSERT(thermal, "thermal BC specification not understood", assert::User_error)
    bc = std::make_shared<No_slip>(thermal, heat_rat, _solver().viscosity_model(), _solver().turbulence_model(),
                                   _vard("heat_flux_coercion"));
  } else if (name == "expression") {
    HEXED_ASSERT(_inter.variables->lookup<std::string>("bc_state"),
                 "To use the `expression` BC type, you must define `bc_state` as a string.",
                 assert::User_error)
    std::string flux;
    if (_inter.variables->lookup<std::string>("bc_flux")) {
      flux = _vars("bc_flux");
    } else {
      for (int i_var = 0; i_var < _vari("n_var"); ++i_var) {
        flux += format_str("ghost_flux%i = flux%i\n", i_var, i_var);
      }
    }
    bc = std::make_shared<Expression_bc>(_inter, _vars("bc_state"), flux);
  } else HEXED_THROW(format_str(1000, "unrecognized boundary condition type `%s`", name.c_str()), assert::User_error);
  return bc;
}

std::string Case::_iteration_suffix() {
  return format_str(100, "iter%.*i", _vari("iter_width"), _vari("iteration"));
}

void force_symlink(const std::filesystem::path& target, const std::filesystem::path& link) {
  if (std::filesystem::exists(link)) std::filesystem::remove(link);
  std::filesystem::create_symlink(target, link);
}

std::vector<std::shared_ptr<Flow_bc>> Case::_make_extremal_bcs() {
  std::vector<std::shared_ptr<Flow_bc>> bcs;
  for (int i_dim = 0; i_dim < _vari("n_dim"); ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      bcs.push_back(_make_bc(_vars(format_str(50, "extremal_bc%i%i", i_dim, sign))));
    }
  }
  return bcs;
}

std::shared_ptr<Surface_geom> Case::_make_geom() {
  int nd = _vari("n_dim");
  Int n_div_min = math::pow(Int(2), _vari("min_geom_subdiv_levels"));
  Int n_div_max = math::pow(Int(2), _vari("max_geom_subdiv_levels"));
  std::vector<Surface_geom*> geoms;
  for (int i_geom = 0;; ++i_geom) {
    auto geom = _inter.variables->lookup<std::string>("geom" + std::to_string(i_geom));
    if (!geom) break;
    HEXED_ASSERT(std::filesystem::exists(geom.value()),
                 format_str(1000, "geometry file `%s` not found", geom->c_str()), assert::User_error)
    Task_message tm(printers::info, "  reading geometry file `" + geom.value() + "`");
    std::string ext = file_extension(geom.value());
    std::string without_ext(geom->begin(), geom->end() - ext.size() - 1);
    for (char& c : ext) c = tolower(c);
    if (ext == "csv") {
      HEXED_ASSERT(nd == 2, "3D geometry in CSV format is not supported", assert::User_error);
      Mat<dyn, dyn> data = read_csv(*geom).transpose();
      HEXED_ASSERT(data.rows() >= nd, "CSV geometry file must have at least n_dim columns", assert::User_error);
      Array<double> data_arr({data.cols(), 3});
      data_arr = 0;
      for (int row = 0; row < data.cols(); ++row) data_arr(row)(0, 2).vector() = data(all, row);
      Tree_curve_geom* geom = new Tree_curve_geom(data_arr.copy(), 4);
      geoms.emplace_back(geom);
    } else if (ext == "igs" || ext == "iges") {
      if (nd == 3) {
        auto ptr = std::make_unique<brep::Geom_3d>(geom.value(), n_div_min, n_div_max, _vard("coincidence_tol_bbox"),
                                                   _vard("coincidence_tol_abs"), _vard("coincidence_tol_subdiv"),
                                                   _vard("tangency_tol_angle"), _vard("tangency_tol_subdiv"));
        if (_vari("vis_geom")) {
          Mat<3, 2> bounds;
          for (int i_dim = 0; i_dim < 3; ++i_dim) {
            for (int sign = 0; sign < 2; ++sign) {
              bounds(i_dim, sign) = _vard("geom_vis_bound" + std::to_string(i_dim) + std::to_string(sign));
            }
          }
          ptr->visualize("default", _vars("working_dir") + without_ext,
                         _vari("geom_vis_subdivisions"), _vari("vis_geom_distance"), bounds);
        }
        geoms.emplace_back(ptr.release());
      } else if (nd == 2) {
        auto ptr = std::make_unique<brep::Geom_2d>(geom.value(), n_div_max);
        if (_vari("vis_geom")) {
          ptr->visualize("default", _vars("working_dir") + without_ext, _vari("geom_vis_subdivisions"));
        }
        geoms.emplace_back(ptr.release());
      } else HEXED_THROW("BRep geometry must be 2 or 3D", assert::User_error)
    } else {
      HEXED_ASSERT(false, format_str(1000, "file extension `%s` not recognized", ext.c_str()), assert::User_error);
    }
  }
  return geoms.empty() ? std::shared_ptr<Surface_geom>() : std::make_shared<Compound_geom>(geoms);
}

std::string Case::_assignment(std::string var_name) {
  std::string statement = var_name + " = ";
  if      (_inter.variables->lookup<        int>(var_name)) statement += std::to_string(_vari(var_name));
  else if (_inter.variables->lookup<     double>(var_name)) statement += format_str(100, "%.20e", _vard(var_name));
  else if (_inter.variables->lookup<std::string>(var_name)) statement += "{" + _vars(var_name) + "}";
  return statement;
}

void Case::_visualize(std::string suffix) {
  Task_message(printers::info, "visualizing");
  std::string wd = _vars("working_dir");
  int n_sample = _vari("vis_n_sample");
  std::vector<std::string> vis_objects {"surface", "field"};
  for (int i_contour = 0; ; ++i_contour) {
    std::string name = "contour" + std::to_string(i_contour);
    if (_inter.variables->lookup<std::string>(name)) vis_objects.push_back(name);
    else break;
  }
  for (std::string v : vis_objects) if (_vari("vis_" + strip_trailing_digits(v))) {
    for (std::string format : {"xdmf", "tecplot", "csv"}) if (_vari("vis_" + format)) {
      std::string vis_expr = _vars("vis_" + strip_trailing_digits(v) + "_vars");
      for (bool edges : {false, true}) {
        std::string name = v;
        if (edges) name = name + "_edges";
        if (!edges || (_vari("vis_edges") && _vari("n_dim") + (v == "field") > 2)) {
          std::string file_name = wd + name + suffix;
          if (v == "surface") {
            _solver().visualize_surface(format, file_name, _solver().mesh().surface_bc_sn(), vis_expr, n_sample, edges);
          } else if (v == "field") {
            _solver().visualize_field(format, file_name, vis_expr, n_sample, edges);
            if (_vari("vis_lts_constraints")) {
              _solver().vis_lts_constraints(format, wd + "lts_constraints" + suffix, n_sample);
            }
          } else if (!edges) { // vis_type == contour0, contour1, etc
            std::string contour_expr = _vars("vis_contour_vars") + v + "_var = " + _vars(v) + ";";
            auto tol = _inter.variables->lookup<double>(name + "_tol");
            double const_tol = (tol ? *tol : 1e-10)*_vard("general_tolerance");
            _solver().visualize_contour(format, file_name, _vars(v), _vars("vis_contour_vars"), const_tol, n_sample);
          }
          if (format == "xdmf") {
            std::string latest = wd + name + "_latest1.xmf";
            if (std::filesystem::exists(latest)) {
              std::filesystem::copy_file(latest, wd + name + "_latest0.xmf",
                                         std::filesystem::copy_options::overwrite_existing);
            }
            if (std::filesystem::exists(file_name + ".xmf")) {
              std::filesystem::copy_file(file_name + ".xmf", latest, std::filesystem::copy_options::overwrite_existing);
            }
          }
        }
      }
    }
  }
}

Transport_model Case::_transport_model(std::string name) {
  auto sub = _inter.make_sub();
  Struct_expr model(_vars(name + "_model"));
  model.eval(sub);
  if (model.names.empty()) {
    return inviscid;
  } else if (sub.variables->exists("offset")) {
    return Transport_model::sutherland(sub.variables->get<double>("ref_value"),
                                       sub.variables->get<double>("ref_temperature"),
                                       sub.variables->get<double>("offset"));
  } else if (sub.variables->exists("const_value")) {
    return Transport_model::constant(sub.variables->get<double>("const_value"));
  } else HEXED_THROW(format_str(200, "invalid transport model specification `{%s}` for %s",
                                model, name.c_str()), assert::User_error) throw;
}

std::function<bool(Element&, int)> Case::_ref_crit(std::string name) {
  return [this, name](Element& elem, int i_dim)->bool {
    auto sub = _inter.make_sub();
    vis_variables::element(*sub.variables, elem);
    sub.variables->assign("i_dim", i_dim);
    sub.exec("return = $" + name);
    return sub.variables->get<int>("return");
  };
}

void Case::_update_monitors() {
  Int iter = _vari("iteration");
  std::ofstream status_data ((_vars("working_dir") + "status_data.txt").c_str());
  auto sub = _inter.make_sub();
  sub.exec(_vars("monitor_vars"));
  _log_residual_hist.add_sample(iter, std::log(_vard("normalized_residual")));
  bool allow_ref = iter > _vari("refine_start_iter");
  if (_vari("automate_adapt_schedule")) {
    bool sufficient_drop = _vard("normalized_residual") < _vard("next_refine_residual");
    bool stagnated = _log_residual_hist.converged({_vard("residual_stagnation_tol")}, {huge});
    bool converged = _vari("monitors_converged");
    allow_ref = allow_ref && (sufficient_drop || stagnated || converged);
  }
  _inter.variables->assign<int>("allow_refinement", allow_ref);
  for (Int i_monitor = -1; i_monitor < Int(_monitor_vars.size()); ++i_monitor) {
    auto& monitor = i_monitor >= 0 ? _monitors[i_monitor] : _log_residual_hist;
    std::string var_name = "log_normalized_residual";
    if (i_monitor >= 0) {
      var_name = _monitor_vars[i_monitor];
      double val = sub.variables->get<double>(var_name);
      monitor.add_sample(iter, val);
    }
    auto assign = [&](std::string suffix, double value) {
      std::string name = var_name + suffix;
      _inter.variables->assign(name, value);
      status_data << name << ": " << to_string(value) << "\n";
    };
    assign("_smoothed", monitor.smoothed());
    assign("_trend", monitor.trend());
    assign("_curvature", monitor.curvature());
    assign("_noise", monitor.noise());
    assign("_noise_trend", monitor.noise_trend());
    assign("_noise_curvature", monitor.noise_curvature());
  }
}

Case::Case(std::string input_script)
: _start_time{std::chrono::system_clock::to_time_t(std::chrono::system_clock::now())}
, _log_residual_hist{.4}
{
  _inter.variables->assign("input_script", input_script);
  _inter.variables->assign("version_major", config::version_major);
  _inter.variables->assign("version_minor", config::version_minor);
  _inter.variables->assign("version_patch", config::version_patch);
  _inter.variables->assign<std::string>("commit", config::commit);
  _inter.variables->assign("vis_default_format", Visualizer::default_format);
  _inter.variables->assign<int>("start_time", _start_time);

  // create custom Heisenberg variables

  _inter.variables->create("setup_output", new Namespace::Heisenberg<std::string>([this]() {
    _output_file.reset(new std::ofstream(_vars("working_dir") + "output.txt"));
    for (auto* printer : {&printers::info, &printers::warn, &printers::error}) {
      printer->printers.emplace_back(std::make_shared<Stream_printer>(*_output_file));
    }
    char utc [100];
    std::strftime(utc, 100, "%Y-%m-%d %H:%M:%S", std::gmtime(&_start_time));
    printers::info(format_str(1000, "Commencing simulation with Hexed version %i.%i.%i (commit %s) at %s UTC (%i Unix Time).\n",
                              config::version_major, config::version_minor, config::version_patch, config::commit.c_str(), utc, _start_time));
    // setup storage parameters
    int n_dim = _inter.variables->get<int, assert::User_error>("n_dim", "User must define `n_dim`.");
    HEXED_ASSERT((n_dim > 0) && (n_dim <= 3), "`n_dim` must be an integer in [1, 3]", assert::User_error);
    HEXED_ASSERT(_inter.variables->exists_recursive("row_size"), "no row size ??");
    int row_size = _inter.variables->get<int, assert::User_error>(
      "row_size",
      "`row_size` is not defined as an integer (did you define it as a different type?)"
    );
    HEXED_ASSERT(row_size >= 2 && row_size <= config::max_row_size,
                 str_cat("`row_size` must be between 2 and ", config::max_row_size), assert::User_error)
    _inter.variables->assign("n_var", n_dim + 2 + 2*(_vars("turbulence_model") == "k-omega"));
    _inter.variables->get<double, assert::User_error>("reference_length", "Must specify `reference_length`.");
    if (n_dim == 3) {
      _inter.variables->get<double, assert::User_error>("reference_area", "Must specify `reference_area`.");
    } else {
      HEXED_ASSERT(!_inter.variables->exists_recursive("reference_area"), "Cannot specify `reference_area` in 2D",
                   assert::User_error)
    }
    return "";
  }));

  _inter.variables->create("compute_freestream", new Namespace::Heisenberg<std::string>([this]() {
    int n_dim = _vari("n_dim");
    int n_var = _vari("n_var");
    Mat<> freestream(n_var);
    if (_inter.variables->lookup<double>("freestream0")) {
      freestream = _get_vector("freestream", n_dim + 2);
    } else {
      if (_inter.variables->lookup<double>("altitude")) {
        HEXED_ASSERT(!_inter.variables->lookup<double>("freestream_temperature"),
                     "cannot specify both altitude and temperature (consider `temperature_offset`)",
                     assert::User_error)
        auto dens_pres = standard_atmosphere(_vard("altitude"), _vard("temperature_offset"));
        _inter.variables->assign<double>("freestream_density", dens_pres[0]);
        _inter.variables->assign<double>("freestream_pressure", dens_pres[1]);
      }
      HEXED_ASSERT(  _inter.variables->lookup<double>("freestream_density").has_value()
                   + _inter.variables->lookup<double>("freestream_pressure").has_value()
                   + _inter.variables->lookup<double>("freestream_temperature").has_value() == 2,
                   "exactly two of freestream density, pressure, and temperature must be specified",
                   assert::User_error)
      if (_inter.variables->lookup<double>("freestream_density")) {
        freestream(n_dim) = _vard("freestream_density");
        if (_inter.variables->lookup<double>("freestream_pressure")) {
          _inter.variables->assign<double>("freestream_temperature",
            _vard("freestream_pressure")/(constants::specific_gas_air*_vard("freestream_density")));
        } else {
          _inter.variables->assign<double>("freestream_pressure",
            _vard("freestream_density")*constants::specific_gas_air*_vard("freestream_temperature"));
        }
      } else {
        _inter.variables->assign<double>("freestream_density",
          _vard("freestream_pressure")/(constants::specific_gas_air*_vard("freestream_temperature")));
      }
      HEXED_ASSERT(  _inter.variables->lookup<double>("freestream_velocity0").has_value()
                   + _inter.variables->lookup<double>("freestream_speed").has_value()
                   + _inter.variables->lookup<double>("freestream_mach").has_value() == 1,
                   "exactly one of velocity, speed, and Mach number must be specified", assert::User_error);
      Mat<> veloc;
      Mat<> full_direction = Mat<>::Zero(3);
      auto direction = full_direction(Eigen::seqN(0, n_dim));
      if (_inter.variables->lookup<double>("freestream_velocity0")) {
        veloc = _get_vector("freestream_velocity", n_dim);
        direction = veloc.normalized();
      } else {
        double fss = std::sqrt(heat_rat*constants::specific_gas_air*_vard("freestream_temperature"));
        _inter.variables->assign<double>("freestream_sound_speed", fss);
        if (_inter.variables->lookup<double>("freestream_speed")) {
          _inter.variables->assign<double>("freestream_mach", _vard("freestream_speed")/_vard("freestream_sound_speed"));
        } else {
          _inter.variables->assign<double>("freestream_speed", _vard("freestream_mach")*_vard("freestream_sound_speed"));
        }
        if (_inter.variables->lookup<double>("freestream_direction0")) {
          direction = _get_vector("freestream_direction", n_dim).normalized();
        } else {
          direction.setUnit(n_dim, 0);
          if (n_dim == 2) {
            direction = Eigen::Rotation2D<double>(_vard("attack"))*direction;
          }
          if (n_dim == 3) {
            direction = Eigen::AngleAxis<double>(-_vard("attack"  ), Eigen::Vector3d::Unit(1))*direction;
            direction = Eigen::AngleAxis<double>( _vard("sideslip"), Eigen::Vector3d::Unit(2))*direction;
          }
        }
        veloc = _vard("freestream_speed")*direction;
        _set_vector("freestream_velocity", veloc);
        for (int i_dim = n_dim; i_dim < 3; ++i_dim) {
          _inter.variables->assign("freestream_velocity" + std::to_string(i_dim), 0.);
        }
      }
      _set_vector("freestream_direction", full_direction);
      double density = _vard("freestream_density");
      double ener = _vard("freestream_pressure")/(heat_rat - 1) + .5*density*veloc.squaredNorm();
      _inter.variables->assign("freestream_energy", ener);
      double dyn_visc = _transport_model("viscosity").coefficient(std::sqrt(_vard("freestream_temperature")));
      _inter.variables->assign("freestream_dynamic_viscosity", dyn_visc);
      double therm_cond = _transport_model("conductivity").coefficient(std::sqrt(_vard("freestream_temperature")));
      _inter.variables->assign("freestream_thermal_conductivity", therm_cond);
      freestream(Eigen::seqN(0, n_dim)) = density*veloc;
      freestream(n_dim) = density;
      freestream(n_dim + 1) = ener;
      if (_transport_model("viscosity").is_viscous) {
        double rpl = density*_vard("freestream_speed")/dyn_visc;
        _inter.variables->assign("reynolds_per_length", rpl);
        double min_shock_width = _vard("min_wall_shock_width_reynolds")/rpl;
        _inter.variables->assign("wall_shock_width", std::max(_vard("wall_shock_width"), min_shock_width));
      } else {
        _inter.variables->assign("reynolds_per_length", huge);
        _inter.variables->assign("reynolds", huge);
      }
      double stke = 0.;
      double turb_diss = 0.;
      double turb_len = 0.;
      if (_vars("turbulence_model") == "k-omega") {
        double speed = _vard("freestream_speed");
        stke = 1.5*math::pow(_vard("freestream_turbulence_intensity")*speed, 2);
        freestream(n_dim + 2) = _vard("freestream_density")*stke;
        turb_len = _vard("freestream_turbulence_length");
        if (turb_len < 0) turb_len = _vard("reference_length");
        turb_diss = 5.*speed/turb_len;
        freestream(n_dim + 3) = _vard("freestream_density")*std::log(turb_diss);
      }
      _inter.variables->assign("freestream_specific_turbulent_kinetic_energy", stke);
      _inter.variables->assign("freestream_turbulence_length", turb_len);
      _inter.variables->assign("freestream_specific_turbulent_dissipation", turb_diss);
      _set_vector("freestream", freestream);
    }
    return "";
  }));

  /*"
   * Creates a `hexed::Solver` object.
   * Mesh contains a single element and the flow is uninitialized.
  "*/
  _inter.variables->create("create_solver", new Namespace::Heisenberg<std::string>([this]() {
    int n_dim = _vari("n_dim");
    double domain_size = _vard("domain_size");
    if (domain_size <= 0) domain_size = _vard("domain_size_ratio")*_vard("reference_length");
    // evaluate dimensions
    Mat<dyn, dyn> mesh_extremes(n_dim, 2);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      for (int sign = 0; sign < 2; ++sign) {
        mesh_extremes(i_dim, sign) = _vard(str_cat("domain_center", i_dim)) + math::sign(sign)*domain_size/2;
      }
    }
    // construct molecular transport models
    std::vector<std::string> transport_phenomena {"viscosity", "conductivity"};
    std::vector<Transport_model> transport_models;
    for (std::string name : transport_phenomena) transport_models.push_back(_transport_model(name));
    Turbulence_model turb_model;
    std::string turb = _vars("turbulence_model");
    if      (turb == "") turb_model = laminar;
    else if (turb == "k-omega") turb_model = k_omega;
    else HEXED_THROW("unrecognized turbulence model `{" + turb + "}`", assert::User_error);
    // setup actual solver
    bool steady = _vari("steady");
    bool implicit = _vari("implicit");
    std::string ts_str = to_lower(_vars("implicit_scheme"));
    Time_scheme ts;
    if (steady) ts = explicit_steady;
    else if (!implicit) ts = explicit_unsteady;
    else if (ts_str == "backward euler") ts = backward_euler;
    else if (ts_str == "crank-nicolson") ts = crank_nicolson;
    else if (ts_str == "dirk2") ts = dirk2;
    else {
      HEXED_THROW("`" + ts_str + "` is not a supported time integration scheme.") throw;
    }
    _solver_ptr.reset(new Solver(n_dim, _vari("row_size"), domain_size, ts, transport_models[0],
                                 transport_models[1], turb_model, _inter.variables));
    _solver().mesh().add_tree(_make_extremal_bcs(), mesh_extremes(all, 0));
    _solver().set_fix_nonphysical(_vari("fix_nonphysical"));
    return "";
  }));

  /*"
   * Creates the mesh, performs initial refinement, and incorporates the surface geometry.
   * You can define 0 or more string variables named `geom0`, `geom1`, `geom2`, ... containing paths to geometry definition files.
   * Paths can be absolute or relative to execution directory (not \ref working_dir).
   * Geometry variable numbers must start at 0 and be consecutive.
   * E.g., if you define `geom0 = leading.txt`, `geom1 = main.txt`, `geom3 = trailing.txt`, then `geom3` will be ignored.
   * If these geometries split the domain into disjoint regions,
   * then the region containing the \ref flood_fill_startN "flood fill start" shall be meshed.
   * Geometry file format is inferred from the file extension, which is case-insensitive.
   * The following extension/format combinations are supported:
   * - `.csv`: Comma Separated Value format (2D only).
   *   A comma-delimited table with at least 2 columns, which are interpreted as \f$ x_0 \f$ and \f$ x_1 \f$ coordinates of the nodes of a polygonal curve,
   *   in \ref units "standard units".
   *   Curve is treated as open---if you are trying to model a closed shape, simply make the first point (approximately) equal to the last.
   * - `.igs`, `.iges`: [IGES](https://en.wikipedia.org/wiki/IGES) CAD
   *   [format](https://filemonger.com/specs/igs/devdept.com/version6.pdf) (2D or 3D).
   *   For 2D simulations, the model curves/edges will be extracted. All curves must (approximately) lie in the \f$ (x_0, x_1) \f$ plane,
   *   and any deviation from said plane will be a source of numerical error.
   *   For 3D simulations, the model surfaces will be used.
   *   Units are read from the file and converted.
   *   Files must be in ASCII format (there is also a binary format, although this is less common).
   *   Both GNU/Linux and Windows line endings are accepted.
  "*/
  _inter.variables->create("mesh", new Namespace::Heisenberg<std::string>([this]() {
    printers::info("Meshing...\n");
    auto compute_bbox = [&]() {
      _solver().bounds_surface("position0 = pos0; position1 = pos1; position2 = pos2;", 2*_vari("n_dim"), 20);
      double geom_len = 0;
      for (int i_dim = 0; i_dim < _vari("n_dim"); ++i_dim) {
        std::vector<std::string> minmax {"min", "max"};
        for (int sign : {0, 1}) {
          _inter.variables->assign("geom_bbox" + to_string(i_dim) + to_string(sign),
                                   _vard(minmax[sign] + "_surface_position" + to_string(i_dim)));
        }
        double dim_len =   _vard("max_surface_position" + to_string(i_dim))
                         - _vard("min_surface_position" + to_string(i_dim));
        geom_len = std::max(geom_len, dim_len);
      }
    };
    auto refine_isotropic = [&](std::string short_name, std::string long_name, bool bbox, bool newline) {
      std::vector<std::string> crit_names {"_refine_if", "_unrefine_if"};
      std::vector<std::function<bool(Element&)>> crits;
      for (std::string crit : crit_names) {
        crits.emplace_back([this, crit, short_name](Element& elem) {
          auto sub = _inter.make_sub();
          vis_variables::element(*sub.variables, elem);
          sub.exec("return = $" + short_name + crit);
          return sub.variables->get<int>("return");
        });
      }
      for (int i_ref = 0, changed = true; i_ref < _vari("max_" + short_name + "_refine_iters") && changed; ++i_ref) {
        printers::info("  " + long_name + " refinement sweep " + to_string(i_ref) + "..." + (newline ? "\n" : " "));
        changed = _solver().mesh().update(crits[0], crits[1]);
        _solver().calc_jacobian();
        if (bbox) compute_bbox();
        printers::info((newline ? "  " : "") + std::string("done. Mesh has ")
                       + to_string(_solver().mesh().n_elements()) + " elements." + (newline ? "\n  " : " "));
        _inter.variables->assign("flow_time", double(i_ref));
        _visualize("_" + short_name + "_ref_sweep" + to_string(i_ref));
      }
    };
    refine_isotropic("init", "Initial", false, false);
    auto geom = _make_geom();
    if (geom.use_count()) {
      printers::info("  Fitting geometry...\n");
      _has_geom = true;
      _solver().mesh().set_surface(geom, _make_bc(_vars("surface_bc")),
                                   _get_vector("flood_fill_start", _vari("n_dim")));
      _solver().calc_jacobian();
      compute_bbox();
        printers::info("  done. Mesh has " + to_string(_solver().mesh().n_elements()) + " elements.\n  ");
      _inter.variables->assign("flow_time", 0.);
      _visualize("_init_geometry_fit");
    }
    refine_isotropic("geom", "Geometry", true, true);
    _inter.variables->assign("flow_time", 0.);
    for (int i_split = 0; i_split < _vari("init_layer_splits"); ++i_split) _inter.make_sub().exec("split_layers");
    _solver().init_wall_dist();
    _solver().update_wall_dist(_vari("wall_dist_iters_initial"));
    for (int i_ref = 0; i_ref < _vari("max_final_refine_iters"); ++i_ref) {
      printers::info("  Final refinement sweep " + to_string(i_ref) + "... ");
      std::vector<std::string> crit_names {"_refine_if", "_unrefine_if"};
      std::vector<std::function<bool(Element&, int)>> crits;
      auto result = _solver().mesh().plan_adaptation(_ref_crit("final_refine_if"),
                                                     _ref_crit("final_unrefine_if"), true);
      if (result.changed) _solver().mesh().execute_adaptation();
      _solver().calc_jacobian();
      _solver().update_wall_dist(_vari("wall_dist_iters_update"));
      printers::info("done. Mesh has " + to_string(_solver().mesh().n_elements()) + " elements. ("
                     + to_string(result.n_refine) + " new " + to_string(result.n_coarsen) + " lost)");
      _inter.variables->assign("flow_time", double(i_ref));
      _visualize("_final_ref_sweep" + to_string(i_ref));
      if (!result.changed) break;
    }
    _inter.variables->assign("hexed_mesh_init", 1);
    printers::info("  geometry bounding box: \n");
    for (int i_dim = 0; i_dim < _vari("n_dim"); ++i_dim) {
      printers::info("   ");
      for (int sign : {0, 1}) {
        printers::info(" " + to_string(_vard("geom_bbox" + to_string(i_dim) + to_string(sign))));
      }
      printers::info("\n");
    }
    printers::info("Meshing complete with " + to_string(_solver().mesh().n_elements()) + " elements.\n", true);
    _solver().update_preti_iters();
    _solver().print_preti_iters();
    return "";
  }));

  _inter.variables->create("refine", new Namespace::Heisenberg<int>([this]() {
    std::vector<std::string> crit_code;
    crit_code.push_back("return = " + _vars("refine_if"));
    crit_code.push_back("return = " + _vars("unrefine_if"));
    std::vector<std::function<bool(Element&)>> crits;
    for (std::string code : crit_code) {
      crits.emplace_back([this, code](Element& elem) {
        auto sub = _inter.make_sub();
        vis_variables::element(*sub.variables, elem);
        sub.exec(code);
        return sub.variables->get<int>("return");
      });
    }
    bool changed = _solver().mesh().update(crits[0], crits[1]);
    _solver().calc_jacobian();
    _visualize("_ref_sweep" + to_string(_vari("i_refinement") + 1));
    return changed;
  }));

  _inter.variables->create("adapt", new Namespace::Heisenberg<std::string>([this]() {
    int iter = _vari("iteration");
    if (iter < _vari("adapt_start_iter") || iter > _vari("adapt_stop_iter")) return "";
    bool allow_ref = _vari("allow_refinement");
    printers::info("Performing uncertainty-based mesh adaptation (");
    if (allow_ref) {
      printers::info("refinement allowed", true);
      _inter.variables->assign("next_refine_residual", _vard("normalized_residual")*_vard("adapt_residual_factor"));
      _inter.variables->assign("last_refine_iter", iter);
    } else {
      printers::info("only coarsening allowed");
    }
    printers::info(")...");
    _solver().compute_spectral_uncertainty();
    double tol_factor = 1;
    Mesh::Adaptation_result result;
    if (allow_ref) {
      while (true) {
        _inter.variables->assign("hexed_tol_factor", tol_factor);
        result = _solver().mesh().plan_adaptation(_ref_crit("adapt_refine_if"), _ref_crit("adapt_unrefine_if"), true);
        if (tol_factor*_vard("general_tolerance")*std::min(_vard("spectral_tol"), _vard("flux_tol")) > 1e3) {
          printers::warn("\n  Only coarsening because excessive refinement could not be avoided "
                         "without excessive tolerance.", true);
          result = _solver().mesh().plan_adaptation([](Element&, int){return false;},
                                                    _ref_crit("adapt_unrefine_if"), true);
          break;
        }
        if (_solver().mesh().n_elements() + result.n_refine - result.n_coarsen < _vard("max_n_elements")) {
          break;
        } else {
          tol_factor *= 1.2;
          printers::warn("\n  Temporarily increasing refinement tolerance to "
                         + to_string(tol_factor*_vard("general_tolerance"))
                         + " to satisfy refinement constraints.", true);
        }
      }
    } else {
      result = _solver().mesh().plan_adaptation([](Element&, int){return false;}, _ref_crit("adapt_unrefine_if"), true);
    }
    _inter.variables->assign("hexed_tol_factor", tol_factor);
    if (result.changed) {
      _solver().mesh().execute_adaptation();
    }
    _solver().calc_jacobian();
    _solver().update_wall_dist(_vari("wall_dist_iters_update"));
    _inter.variables->assign<int>("adapt_changed", result.changed);
    _solver().compute_residual(false);
    printers::info(" done. Mesh now has " + to_string(_solver().mesh().n_elements()) + " elements.\n");
    _solver().update_preti_iters();
    _solver().print_preti_iters();
    return "";
  }));

  _inter.variables->create("adapt_shock", new Namespace::Heisenberg<std::string>([this]() {
    printers::info("Performing quick-response adaptation (artificial viscosity-based)...\n");
    for (Int sweep = 0; sweep < _vari("shock_refine_iters"); ++sweep) {
      printers::info("  coarsening sweep " + to_string(sweep) + ":");
      auto result = _solver().mesh().plan_adaptation([](Element&, int){return false;},
                                                     _ref_crit("shock_unrefine_if"), false);
      if (result.changed) _solver().mesh().execute_adaptation();
      printers::info(" " + to_string(_solver().mesh().n_elements()) + " elements\n");
      if (!result.changed) break;
    }
    for (Int sweep = 0; sweep < _vari("shock_refine_iters"); ++sweep) {
      printers::info("  refinement sweep " + to_string(sweep) + ":");
      auto result = _solver().mesh().plan_adaptation(_ref_crit("shock_refine_if"),
                                                     [](Element&, int){return false;}, false);
      if (result.n_elements > _vard("max_n_elements")) {
        printers::error(" aborting refinement to avoid exceeding maximum number of elements!", true);
        break;
      }
      if (result.changed) _solver().mesh().execute_adaptation();
      printers::info(" " + to_string(_solver().mesh().n_elements()) + " elements\n");
      if (!result.changed) break;
    }
    _solver().calc_jacobian();
    _solver().update_preti_iters();
    _solver().print_preti_iters();
    _solver().compute_residual(false);
    printers::info("done\n");
    return "";
  }));

  /*" Initializes the state to \ref init_cond "*/
  _inter.variables->create("init_state", new Namespace::Heisenberg<std::string>([this]() {
    _solver().initialize(_vars("init_cond"));
    // implicit setup
    bool implicit = !_vari("steady") && _vari("implicit");
    if (implicit) {
      HEXED_ASSERT(_inter.variables->lookup<double>("time_step"),
                   "unsteady implicit time marching requires you to set `time_step` to a floating-point value.",
                   assert::User_error)
      HEXED_ASSERT(_vard("time_step") >= 0, "`time_step` must be nonnegative.", assert::User_error)
      _inter.variables->assign("flow_time", _vard("flow_time") + _vard("time_step"));
      _inter.variables->assign("hexed_next_flow_time", _vard("flow_time") + _vard("time_step"));
    }
    return "";
  }));

  _inter.variables->create("init_monitors", new Namespace::Heisenberg<std::string>([this]() {
    auto sub = _inter.make_sub();
    sub.subspace();
    sub.exec(_vars("monitor_vars"));
    _monitor_vars = sub.variables->names();
    for (std::string name : _monitor_vars) {
      _monitors.emplace_back(_vard("monitor_window"));
      for (std::string suffix : {"_smoothed", "_noise", "_trend", "_curvature", "_noise_trend", "_noise_curvature"}) {
        _inter.variables->assign_default(name + suffix, 0.);
      }
    }
    return "";
  }));

  /*"
   * Reads the mesh from a Hexed mesh file (which you can create using \ref write_mesh).
   * This can be used in place of \ref mesh, but you will still have to call \ref init_state or \ref read_state to initialize the flow state.
   * The file name of the mesh file given by \ref input_data with `.mesh.h5` appended.
  "*/
  _inter.variables->create("read_mesh", new Namespace::Heisenberg<std::string>([this]() {
    printers::info("Reading mesh...\n");
    auto geom = _make_geom();
    _solver().read_mesh(_vars("input_data"), _make_extremal_bcs(), geom,
                        geom.use_count() ? _make_bc(_vars("surface_bc")) : nullptr);
    _solver().update_preti_iters();
    _solver().print_preti_iters();
    printers::info("done\n");
    return "";
  }));

  /*"
   * Reads the flow state from a file created with the \ref write_state command.
   * The file name is given by \ref input_data with `.state.h5` appended.
   * The state file __must__ have been created for the exact same mesh,
   * meaning that if the state file was not created during the same simulation you're currently running,
   * you need to first use \ref read_mesh to get the mesh.
   * Two meshes that look the same might not actually be the same for a variety of reasons.
   * For example, if elements were refined in a different order, then you cannot use the same state files,
   * even if the mesh is geometrically identical.
   * So, just play it safe and read the mesh first.
   * This command can be used instead of \ref init_state.
   * `read_mesh` and `read_state` together can be used to restart a previous simulation from where it left off.
  "*/
  _inter.variables->create("read_state", new Namespace::Heisenberg<std::string>([this]() {
    Task_message(printers::info, "Reading state...\n");
    _solver().read_state(_vars("input_data"));
    printers::info("done\n");
    return "";
  }));

  /*" Reads HIL variables from a file created with the \ref write_status command. "*/
  _inter.variables->create("read_status", new Namespace::Heisenberg<std::string>([this]() {
    printers::info("Reading status...\n");
    auto sub = _inter.make_sub();
    sub.exec("$read {" + _vars("input_data") + ".status.hil}");
    printers::info("done\n");
    return "";
  }));

  /*"
   * Writes the mesh to a file in the Hexed mesh format,
   * which can then be used again in a future simulation using \ref read_mesh, or, in theory, imported into another program.
   * The file name will be `iterXXXXXXXXXX.mesh.h5` in the working directory, where the `X`s are replaced with the current iteration number.
   * Also creates a symlink `latest.mesh.h5` pointing to this file for convenience.
  "*/
  _inter.variables->create("write_mesh", new Namespace::Heisenberg<std::string>([this]() {
    Task_message(printers::info, "writing mesh");
    std::string file_name = _vars("working_dir") + _iteration_suffix();
    _solver().mesh().write(file_name);
    for (std::string suffix : {"mesh", "tree"}) {
      force_symlink(_iteration_suffix() + "." + suffix + ".h5", _vars("working_dir") + "latest." + suffix + ".h5");
    }
    return "";
  }));

  /*"
   * Writes the flow state to a file which can then be used in a future simulation using \ref read_state.
   * The file name will be `iterXXXXXXXXXX.state.h5` in the working directory, where the `X`s are replaced with the current iteration number.
   * Also creates a symlink `latest.state.h5` pointing to this file for convenience.
  "*/
  _inter.variables->create("write_state", new Namespace::Heisenberg<std::string>([this]() {
    Task_message(printers::info, "writing state");
    std::string file_name = _vars("working_dir") + _iteration_suffix();
    _solver().write_state(file_name);
    force_symlink(_iteration_suffix() + ".state.h5", _vars("working_dir") + "latest.state.h5");
    return "";
  }));

  /*"
   * Writes the flow state to a file which can then be used in a future simulation using \ref read_state.
   * The file name will be `iterXXXXXXXXXX.state.h5` in the working directory, where the `X`s are replaced with the current iteration number.
   * Also creates a symlink `latest.state.h5` pointing to this file for convenience.
  "*/
  _inter.variables->create("write_status", new Namespace::Heisenberg<std::string>([this]() {
    Task_message(printers::info, "writing status");
    std::ofstream status_file(_vars("working_dir") + _iteration_suffix() + ".status.hil");
    std::vector<std::string> no_write {"working_dir", "input_data"};
    for (std::string name : _inter.variables->names()) {
      if (name.substr(0, 6) != "hexed_" &&
          std::none_of(no_write.begin(), no_write.end(), [name](std::string nw){return name == nw;})) {
        status_file << _assignment(name) + "\n";
      }
    }
    status_file.close();
    force_symlink(_iteration_suffix() + ".status.hil", _vars("working_dir") + "latest.status.hil");
    return "";
  }));

  /*" Writes visualization files. "*/
  _inter.variables->create("visualize", new Namespace::Heisenberg<std::string>([this]() {
    _visualize("_" + _iteration_suffix());
    printers::info("Current wall clock time: " + to_string(_vard("wall_time")) + "\n");
    return "";
  }));

  _inter.variables->create("write_skews", new Namespace::Heisenberg<std::string>([this]() {
    Csv csv(_vars("working_dir") + "skews", 1);
    Array<double> skews(_solver().skews());
    Array<double> reshaped({skews.size(), 1}, skews.data());
    csv.write(reshaped);
    return "";
  }));

  _inter.variables->create<std::string>("header", new Namespace::Heisenberg<std::string>([this]() {
    std::string header = "";
    auto sub = _inter.make_sub();
    sub.subspace();
    std::string print_expr = _vars("print_vars");
    sub.exec(print_expr);
    _print_vars = sub.variables->names();
    std::sort(_print_vars.begin(), _print_vars.end(), [print_expr](std::string s0, std::string s1) {
      return print_expr.find(s0) < print_expr.find(s1);
    });
    for (std::string name : _print_vars) {
      int width = std::max<int>(15, name.size());
      header += format_str(1000, "%*s, ", width, name.c_str());
    }
    header.erase(header.end() - 2, header.end());
    return header;
  }));

  /*"
   * Evaluates the residuals and assigns them to float variables
   * \ref residual_momentum, \ref residual_density, \ref residual_energy.
   * Also updates the spectral uncertainty in each element.
  "*/
  _inter.variables->create("compute_residuals", new Namespace::Heisenberg<std::string>([this]() {
    int nd = _solver().storage_params().n_dim;
    Physical_residual phys_resid;
    auto compute_res = [&](bool unsteady_implicit, std::string prefix) {
      _solver().compute_residual(unsteady_implicit);
      auto res = _solver().integral_field(Pow(phys_resid, 2));
      for (int i_dim = 1; i_dim < nd; ++i_dim) res[0] += res[i_dim];
      for (double& r : res) r = std::sqrt(r);
      _inter.variables->assign(prefix + "residual_momentum", res[0]);
      _inter.variables->assign(prefix + "residual_density", res[nd]);
      _inter.variables->assign(prefix + "residual_energy", res[nd + 1]);
    };
    compute_res(false, "");
    if (!_vari("steady") && _vari("implicit")) compute_res(true, "unsteady_");
    if (_vari("iteration") > 0) _solver().compute_spectral_uncertainty();
    return "";
  }));

  _inter.variables->create("report", new Namespace::Heisenberg<std::string>([this]() {
    std::string report = "";
    auto sub = _inter.make_sub();
    sub.exec(_vars("print_vars"));
    for (std::string name : _print_vars) {
      int width = std::max<int>(15, name.size());
      std::optional<int> vali;
      std::optional<double> vald;
      std::optional<std::string> vals;
      if ((vali = sub.variables->lookup<int>(name))) {
        report += format_str(1000, "%*i, ", width, vali.value());
      } else if ((vald = sub.variables->lookup<double>(name))) {
        report += format_str(1000, "%*.8e, ", width, vald.value());
      } else if ((vals = sub.variables->lookup<std::string>(name))) {
        report += format_str(1000, "%*s, ", width, vals.value());
      }
    }
    report.erase(report.end() - 2, report.end());
    _solver().reset_counters();
    return report;
  }));

  _inter.variables->create<int>("monitors_converged", new Namespace::Heisenberg<int>([this]() {
    bool converged = _monitor_vars.size();
    double trend_tol = _vard("trend_tol")*_vard("general_tolerance");
    double noise_tol = _vard("noise_tol")*_vard("general_tolerance");
    for (unsigned i_monitor = 0; i_monitor < _monitor_vars.size() && converged; ++i_monitor) {
      converged = converged && _monitors[i_monitor].converged({.rel=trend_tol}, {.rel=noise_tol});
    }
    return converged;
  }));

  /*"
   * Executes \ref print_freq solver iterations and returns an empty string.
   * Each "solver iteration" consists of:
   * -# Calling `hexed::Solver::set_art_visc_constant`, if applicable.
   * -# Calling `hexed::Solver::update_art_visc_smoothness`, if applicable.
   *    This will itself update the artificial viscosity advection equations \ref av_advect_iters times
   *    and the diffusion equations \ref av_diff_iters times.
   * -# Calling `hexed::Solver::update`.
   *    This will itself update the flow equations \ref flow_iters times
   *    (where \ref flow_iters defaults to 1 for simulations without shock-capturing).
  "*/
  _inter.variables->create("update", new Namespace::Heisenberg<std::string>([this]() {
    HEXED_ASSERT(_vari("hexed_mesh_init"), "attempt to update flow when mesh has not been created", assert::User_error);
    _inter.variables->assign("total_smear_iters", 0);
    bool avw = _vari("capture_shocks");
    bool avc = _vard("art_visc_constant") > 0;
    bool be = !_vari("steady") && _vari("implicit");
    int iter = _vari(be ? "pseudotime_iteration" : "iteration");
    int print_freq = _vari("print_freq");
    int n = iter ? print_freq - iter%print_freq : 1;
    for (int i = 0; i < n; ++i) {
      ++iter;
      if (_inter.variables->get<int>("elementwise_art_visc")) {
        HEXED_THROW("Elementwise artificial viscosity is no longer supported.", assert::Not_implemented_error)
      } else if (avw) {
        _solver().update_art_visc_smoothness();
      } else if (avc) {
        _solver().set_art_visc_constant(_vard("art_visc_constant"));
      } else if (_solver().using_art_visc()) {
        printers::info("Turning off artificial viscosity.\n", true);
        _solver().set_art_visc_off();
      }
      _solver().update();
    }
    _inter.variables->assign(be ? "pseudotime_iteration" : "iteration", iter);
    _inter.variables->assign("total_iteration", _vari("total_iteration") + n);
    if (!be) _update_monitors();
    return "";
  }));

  _inter.variables->create<std::string>("next_time_stage", new Namespace::Heisenberg<std::string>([this]() {
    int stage = _solver().next_time_stage();
    _inter.variables->assign("pseudotime_iteration", 0);
    _inter.variables->assign("time_stage", stage);
    if (!stage) {
      _inter.variables->assign("iteration", _vari("iteration") + 1);
      _inter.variables->assign("flow_time", _vard("hexed_next_flow_time"));
      _inter.variables->assign("hexed_next_flow_time", _vard("hexed_next_flow_time") + _vard("time_step"));
      _update_monitors();
    }
    return "";
  }));

  _inter.variables->create<int>("n_elements", new Namespace::Heisenberg<int>([this]() {
    return _solver().mesh().n_elements();
  }));
  _inter.variables->create<std::string>("performance_report", new Namespace::Heisenberg<std::string>([this]() {
    return _solver().stopwatch_tree().report() + _solver().mesh().stopwatch_tree().report();
  }));

  _inter.variables->create<std::string>("update_roughness", new Namespace::Heisenberg<std::string>([this]() {
    _solver().update_bound_conds();
    return "";
  }));

  _inter.variables->create<std::string>("bounds_surface", new Namespace::Heisenberg<std::string>([this]() {
    _solver().bounds_surface(_vars("bounds_surface_vars"), _solver().mesh().surface_bc_sn(), _vari("vis_n_sample"));
    return "";
  }));

  /*"
   * Computes integrals of the variables in \ref integrand_field over the domain
   * by the native quadrature rule of the numerical scheme.
   * The integrals of the variables are assigned to float variables named
   * `integral_field_VAR` where `VAR` is the name of the variable.
  "*/
  _inter.variables->create("integrate_field", new Namespace::Heisenberg<std::string>([this]() {
    _solver().integrate_field(_vars("integrand_field"));
    return "";
  }));

  /*"
   * Computes integrals of the variables in \ref integrand_surface over the geometry surface
   * by the native quadrature rule of the numerical scheme.
   * The integrals of the variables are assigned to float variables
   * named `integral_surface_var` where `var` is the name of the variable.
  "*/
  _inter.variables->create("integrate_surface", new Namespace::Heisenberg<std::string>([this]() {
    _solver().integrate_surface(_vars("integrand_surface"), _solver().mesh().surface_bc_sn());
    return "";
  }));

  // load HIL code for the Case _interface
  _inter.exec("$read {hexed.hil}");
  // make a sub-namespace for all new user-created variables which is useful for `write_status`
  std::shared_ptr<Namespace> space = std::make_shared<Namespace>();
  space->supers.push_back(_inter.variables);
  _inter.variables = space;
  // execute input file
  try {
    try {
      _inter.exec(format_str(1000, "$read {%s}", input_script.c_str()));
    } catch (const assert::Numerical_exception& except) {
      printers::error("\nTerminating simulation due to numerical exception.\n", true);
      if (_solver_ptr) _inter.exec("write_mesh; write_state; write_status; visualize; println performance_report;");
      _inter.variables->assign("hexed_failed", 1);
      throw except;
    }
  } catch (const assert::Exception& except) {
    printers::error("\n" + except.name() + ": ", true);
    printers::error(except.message() + "\n");
  }
}

Case::~Case() {
  for (auto* printer : {&printers::info, &printers::warn, &printers::error}) {
    printer->printers.pop_back();
  }
}

}
