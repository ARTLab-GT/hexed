#include <filesystem>
#include <Case.hpp>
#include <Simplex_geom.hpp>
#include <read_csv.hpp>
#include <standard_atmosphere.hpp>
#include <Occt.hpp>
#include <hil_properties.hpp>

namespace hexed
{

const double heat_rat = 1.4;

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

Flow_bc* Case::_make_bc(std::string name)
{
  Mat<> freestream = _get_vector("freestream", _vari("n_dim").value() + 2);
  if (name == "characteristic") return new Riemann_invariants(freestream);
  else if (name == "freestream") return new Freestream(freestream);
  else if (name == "pressure_outflow") return new Pressure_outflow(_vard("pressure").value());
  else if (name == "nonpenetration") return new Nonpenetration;
  else if (name == "no_slip") {
    auto sub = _inter.make_sub();
    sub.exec("$thermal_bc");
    std::shared_ptr<Thermal_bc> thermal;
    if (sub.variables->exists("heat_flux")) { // note not recursive
      thermal = std::make_shared<Prescribed_heat_flux>(sub.variables->lookup<double>("heat_flux").value());
    } else if (sub.variables->exists("emissivity") || sub.variables->exists("conduction")) {
      auto equilibrium = std::make_shared<Thermal_equilibrium>();
      HEXED_ASSERT(sub.variables->exists("conduction") == sub.variables->exists("temperature"), "must specify both surface conduction and temperature or neigher");
      if (sub.variables->exists("emissivity")) equilibrium->emissivity = sub.variables->lookup<double>("emissivity").value();
      if (sub.variables->exists("conduction")) {
        equilibrium->conduction = sub.variables->lookup<double>("conduction").value();
        equilibrium->temperature = sub.variables->lookup<double>("temperature").value();
      }
      thermal = equilibrium;
    } else if (sub.variables->exists("internal_energy")) { // note not recursive
      thermal = std::make_shared<Prescribed_energy>(sub.variables->lookup<double>("internal_energy").value());
    } else if (sub.variables->exists("temperature")) { // note not recursive
      double energy = sub.variables->lookup<double>("temperature").value()*constants::specific_gas_air/(heat_rat - 1.);
      thermal = std::make_shared<Prescribed_energy>(energy);
    }
    HEXED_ASSERT(thermal, "thermal BC specification not understood");
    return new No_slip(thermal, _vard("heat_flux_coercion").value());
  }
  else HEXED_ASSERT(false, format_str(1000, "unrecognized boundary condition type `%s`", name.c_str()));
  return nullptr; // will never happen. just to shut up GCC warning
}

Case::Case(std::string input_file)
{
  _inter.variables->assign("input_file", input_file);
  _inter.variables->assign("version_major", config::version_major);
  _inter.variables->assign("version_minor", config::version_minor);
  _inter.variables->assign("version_patch", config::version_patch);
  _inter.variables->assign<std::string>("commit", config::commit);
  // create custom Heisenberg variables
  _inter.variables->create("create_solver", new Namespace::Heisenberg<int>([this]() {
    // basic IO setup
    _output_file.reset(new std::ofstream(_vars("working_dir").value() + "output.txt"));
    auto printer = std::make_shared<Stream_printer>();
    printer->streams.emplace_back(_output_file.get());
    _inter.printer = printer;
    char utc [100];
    std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::strftime(utc, 100, "%Y-%m-%d %H:%M:%S", std::gmtime(&time));
    printer->print(format_str(1000, "Commencing simulation with Hexed version %i.%i.%i (commit %s) at %s UTC (%i Unix Time).\n",
                              config::version_major, config::version_minor, config::version_patch, config::commit, utc, time));
    // setup actual solver
    auto n_dim = _inter.variables->lookup<int>("n_dim");
    HEXED_ASSERT(n_dim && n_dim.value() > 0 && n_dim.value() <= 3,
                 "`n_dim` must be defined as an integer in [1, 3]");
    auto row_size = _inter.variables->lookup<int>("row_size");
    HEXED_ASSERT(row_size.value() >= 2 && row_size.value() <= config::max_row_size,
                 format_str(300, "`row_size` must be between 2 and %i", config::max_row_size));
    // compute freestream
    Mat<> freestream(*n_dim + 2);
    if (_vard("freestream0")) freestream = _get_vector("freestream", *n_dim + 2);
    else {
      if (_vard("altitude")) {
        HEXED_ASSERT(!_vard("temperature"), "cannot specify both altitude and temperature (consider `temperature_offset`)");
        auto dens_pres = standard_atmosphere(_vard("altitude").value(), _vard("temperature_offset").value());
        _inter.variables->assign<double>("density", dens_pres[0]);
        _inter.variables->assign<double>("pressure", dens_pres[1]);
      }
      HEXED_ASSERT(_vard("density").has_value() + _vard("pressure").has_value() + _vard("temperature").has_value() == 2,
                   "exactly two of density, pressure, and temperature must be specified");
      if (_vard("density")) {
        freestream(*n_dim) = *_vard("density");
        if (_vard("pressure")) _inter.variables->assign<double>("temperature", *_vard("pressure")/(constants::specific_gas_air**_vard("density")));
        else _inter.variables->assign<double>("pressure", *_vard("density")*constants::specific_gas_air**_vard("temperature"));
      } else _inter.variables->assign<double>("density", *_vard("pressure")/(constants::specific_gas_air**_vard("temperature")));
      HEXED_ASSERT(_vard("velocity0").has_value() + _vard("speed").has_value() + _vard("mach").has_value() == 1,
                   "exactly one of velocity, speed, and Mach number must be specified");
      Mat<> veloc;
      Mat<> full_direction = Mat<>::Zero(3);
      auto direction = full_direction(Eigen::seqN(0, *n_dim));
      if (_vard("velocity0")) {
        veloc = _get_vector("velocity", *n_dim);
        direction = veloc.normalized();
      } else {
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
          _inter.variables->assign<double>("sound", std::sqrt(heat_rat*constants::specific_gas_air**_vard("temperature")));
          if (_vard("speed")) _inter.variables->assign<double>("mach", *_vard("speed")/ *_vard("sound"));
          else _inter.variables->assign<double>("speed", *_vard("mach")**_vard("sound"));
        }
        veloc = *_vard("speed")*direction;
        _set_vector("velocity", veloc);
      }
      _set_vector("direction", full_direction);
      freestream(Eigen::seqN(0, *n_dim)) = *_vard("density")*veloc;
      freestream(*n_dim) = *_vard("density");
      freestream(*n_dim + 1) = *_vard("pressure")/(heat_rat - 1) + .5**_vard("density")*veloc.squaredNorm();
      freestream.conservativeResize(5);
      freestream(Eigen::seqN(*n_dim + 2, 5 - (*n_dim + 2))).setZero();
      _set_vector("freestream", freestream);
    }
    // create solver
    Mat<dyn, dyn> mesh_extremes(*n_dim, 2);
    std::vector<Flow_bc*> bcs;
    for (int i_dim = 0; i_dim < *n_dim; ++i_dim) {
      for (int sign = 0; sign < 2; ++sign) {
        std::string index = format_str(50, "%i%i", i_dim, sign);
        mesh_extremes(i_dim, sign) = _vard("mesh_extreme" + index).value();
        bcs.push_back(_make_bc(_vars("extremal_bc" + index).value()));
      }
    }
    HEXED_ASSERT((mesh_extremes(all, 1) - mesh_extremes(all, 0)).minCoeff() > 0, "all mesh dimensions must be positive!");
    double root_sz = (mesh_extremes(all, 1) - mesh_extremes(all, 0)).maxCoeff();
    std::unique_ptr<Transport_model> visc_model; // make these pointers since assignment operator is deleted
    std::unique_ptr<Transport_model> therm_model;
    if (_vars("transport_model").value() == "inviscid") {
      visc_model.reset(new Transport_model(inviscid));
      therm_model.reset(new Transport_model(inviscid));
    } else if (_vars("transport_model").value() == "sutherland") {
      visc_model.reset(new Transport_model(air_sutherland_dyn_visc));
      therm_model.reset(new Transport_model(air_sutherland_therm_cond));
    } else {
      HEXED_ASSERT(false, "unrecognized transport model specification");
    }
    _solver_ptr.reset(new Solver(*n_dim, *row_size, root_sz, true, *visc_model, *therm_model, _inter.variables, printer));
    _solver().mesh().add_tree(bcs, mesh_extremes(all, 0));
    _solver().set_fix_admissibility(_vari("fix_therm_admis").value());
    return 0;
  }));

  _inter.variables->create<int>("init_refinement", new Namespace::Heisenberg<int>([this]() {
    for (int i = 0; i < _vari("init_ref_level"); ++i) _solver().mesh().update();
    _solver().calc_jacobian();
    return 0;
  }));

  _inter.variables->create<int>("add_geom", new Namespace::Heisenberg<int>([this]() {
    int nd = _vari("n_dim").value();
    std::vector<Surface_geom*> geoms;
    for (int i_geom = 0;; ++i_geom) {
      auto geom = _vars("geom" + std::to_string(i_geom));
      if (!geom) break;
      unsigned dot = geom->rfind('.');
      HEXED_ASSERT(dot < geom->size(), "file name must contain extension to infer format");
      HEXED_ASSERT(std::filesystem::exists(geom.value()), format_str(1000, "geometry file `%s` not found", geom->c_str()));
      std::string case_sensitive(geom->begin() + dot + 1, geom->end());
      std::string ext = case_sensitive;
      for (char& c : ext) c = tolower(c);
      if (ext == "csv") {
        HEXED_ASSERT(nd == 2, "3D geometry in CSV format is not supported");
        auto data = read_csv(*geom);
        HEXED_ASSERT(data.cols() >= nd, "CSV geometry file must have at least n_dim columns");
        geoms.emplace_back(new Simplex_geom<2>(segments(data.transpose())));
      #if HEXED_USE_OCCT
      } else if (ext == "igs" || ext == "iges" || ext == "stp" || ext == "step") {
        auto shape = Occt::read(*geom);
        if (nd == 2) {
          geoms.emplace_back(new Simplex_geom<2>(Occt::segments(shape, _vari("geom_n_segments").value())));
        } else if (nd == 3) {
          auto ptr = new Simplex_geom<3>(Occt::triangles(shape, _vard("max_angle").value(), _vard("max_deflection").value()));
          ptr->visualize(format_str(1000, "%sgeom%i_triangulation", _vars("working_dir").value().c_str(), i_geom));
          geoms.emplace_back(ptr);
        }
      } else if (ext == "stl") {
        HEXED_ASSERT(nd == 3, "STL format is only supported for 3D");
        geoms.emplace_back(new Simplex_geom<3>(Occt::triangles(Occt::read_stl(geom.value()))));
      #endif
      } else {
        HEXED_ASSERT(false, format_str(1000, "file extension `%s` not recognized", case_sensitive.c_str()));
      }
    }
    if (!geoms.empty()) {
      _has_geom = true;
      _solver().mesh().set_surface(new Compound_geom(geoms), _make_bc(_vars("surface_bc").value()), _get_vector("flood_fill_start", nd));
      for (int i_smooth = 0; i_smooth < _vari("n_smooth"); ++i_smooth) _solver().mesh().relax(0.7);
      _solver().calc_jacobian();
    }
    return 0;
  }));

  _inter.variables->create<int>("refine", new Namespace::Heisenberg<int>([this]() {
    std::vector<std::string> crit_code;
    crit_code.push_back("return = " + _vars("surface_refine").value());
    crit_code.push_back("return = " + _vars("surface_unrefine").value());
    std::vector<std::function<bool(Element&)>> crits;
    for (std::string code : crit_code) {
      crits.emplace_back([this, code](Element& elem) {
        auto sub = _inter.make_sub();
        hil_properties::element(*sub.variables, elem);
        sub.exec(code);
        return sub.variables->lookup<int>("return").value();
      });
    }
    Jac_inv_det_func jidf;
    _solver().set_resolution_badness(Elem_nonsmooth(jidf));
    _solver().mesh().set_unref_locks(criteria::if_extruded);
    bool changed = _solver().mesh().update(crits[0], crits[1]);
    for (int i_smooth = 0; i_smooth < _vari("n_smooth"); ++i_smooth) _solver().mesh().relax(0.7);
    _solver().calc_jacobian();
    return changed;
  }));

  _inter.variables->create<int>("split_layers", new Namespace::Heisenberg<int>([this]() {
    _solver().mesh().disconnect_boundary(_solver().mesh().surface_bc_sn());
    auto sub = _inter.make_sub();
    std::vector<double> split_points = Struct_expr(_vars("layer_split_points").value()).eval(sub);
    double prev_split = 1.;
    for (double split : split_points) {
      _solver().mesh().extrude(true, split/prev_split, true);
      prev_split = split;
    }
    _solver().mesh().connect_rest(_solver().mesh().surface_bc_sn());
    _solver().calc_jacobian();
    return 0;
  }));

  _inter.variables->create<int>("init_state", new Namespace::Heisenberg<int>([this]() {
    _solver().initialize(Spacetime_expr(Struct_expr(_vars("init_cond").value()), _inter));
    return 0;
  }));

  _inter.variables->create<int>("visualize", new Namespace::Heisenberg<int>([this]() {
    std::string wd = _vars("working_dir").value();
    std::string suffix = format_str(100, "_iter%.*i", _vari("iter_width").value(), _vari("iteration").value());
    int n_sample = _vari("vis_n_sample").value();
    if (_vari("vis_field").value()) {
      Struct_expr vis_vars(_vars("vis_field_vars").value());
      Qpoint_expr func(vis_vars, _inter);
      std::string file_name = wd + "field" + suffix;
      #if HEXED_USE_TECPLOT
      if (_vari("vis_tecplot").value()) _solver().visualize_field_tecplot(func, file_name, n_sample);
      #endif
      #if HEXED_USE_XDMF
      if (_vari("vis_xdmf").value()) {
        _solver().visualize_field_xdmf(func, file_name, n_sample);
        if (_vari("vis_lts_constraints").value()) _solver().vis_lts_constraints(wd + "lts" + suffix, n_sample);
      }
      #endif
    }
    if (_vari("vis_surface").value() && _has_geom) {
      Struct_expr vis_vars(_vars("vis_surface_vars").value());
      Boundary_expr func(vis_vars, _inter);
      std::string file_name = wd + "surface" + suffix;
      int bc_sn = _solver().mesh().surface_bc_sn();
      #if HEXED_USE_TECPLOT
      if (_vari("vis_tecplot").value()) _solver().visualize_surface_tecplot(bc_sn, func, file_name, n_sample);
      #endif
      #if HEXED_USE_XDMF
      if (_vari("vis_xdmf").value()) _solver().visualize_surface_xdmf(bc_sn, func, file_name, n_sample);
      #endif
    }
    return 0;
  }));

  _inter.variables->create<std::string>("header", new Namespace::Heisenberg<std::string>([this]() {
    std::string header = "";
    Struct_expr vars(_vars("print_vars").value());
    for (std::string name : vars.names) {
      header += format_str(1000, "%14s, ", name.c_str());
    }
    header.erase(header.end() - 2, header.end());
    return header;
  }));

  _inter.variables->create<std::string>("compute_residuals", new Namespace::Heisenberg<std::string>([this]() {
    int nd = _solver().storage_params().n_dim;
    Physical_update update;
    auto res = _solver().integral_field(Pow(update, 2));
    for (double& r : res) r /= math::pow(_inter.variables->lookup<double>("time_step").value(), 2);
    double res_mmtm = 0;
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      res_mmtm += res[i_dim];
    }
    _inter.variables->assign("residual_momentum", res_mmtm);
    _inter.variables->assign("residual_density", res[nd]);
    _inter.variables->assign("residual_energy", res[nd + 1]);
    return "";
  }));

  _inter.variables->create<std::string>("report", new Namespace::Heisenberg<std::string>([this]() {
    _inter.variables->lookup<std::string>("compute_residuals");
    std::string report = "";
    Struct_expr vars(_vars("print_vars").value());
    auto sub = _inter.make_sub();
    for (unsigned i_var = 0; i_var < vars.names.size(); ++i_var) {
      int width = std::max<int>(15, vars.names[i_var].size());
      sub.exec(vars.names[i_var] + " = " + vars.exprs[i_var]);
      std::optional<int> vali;
      std::optional<double> vald;
      std::optional<std::string> vals;
      if ((vali = sub.variables->lookup<int>(vars.names[i_var]))) {
        report += format_str(1000, "%*i, ", width, vali.value());
        _inter.variables->assign(vars.names[i_var], vali.value());
      } else if ((vald = sub.variables->lookup<double>(vars.names[i_var]))) {
        report += format_str(1000, "%*.8e, ", width, vald.value());
        _inter.variables->assign(vars.names[i_var], vald.value());
      } else if ((vals = sub.variables->lookup<std::string>(vars.names[i_var]))) {
        report += format_str(1000, "%*s, ", width, vals.value());
        _inter.variables->assign(vars.names[i_var], vals.value());
      }
    }
    report.erase(report.end() - 2, report.end());
    _solver().reset_counters();
    return report;
  }));

  _inter.variables->create<std::string>("update", new Namespace::Heisenberg<std::string>([this]() {
    bool avw = _vard("art_visc_width").value() > 0;
    bool avc = _vard("art_visc_constant").value() > 0;
    int iter = _vari("iteration").value();
    int print_freq = _vari("print_freq").value();
    int n = iter ? print_freq - iter%print_freq : 1;
    for (int i = 0; i < n; ++i) {
      int iter = _vari("iteration").value();
      if (iter > 0) _solver().update_implicit();
      if (avw) {
        _solver().set_art_visc_smoothness(_vard("art_visc_width").value());
      } else if (avc) {
        _solver().set_art_visc_smoothness(_vard("art_visc_constant").value());
      }
      _solver().update();
      _inter.variables->assign("iteration", iter + 1);
    }
    return "";
  }));
  _inter.variables->create<int>("n_elements", new Namespace::Heisenberg<int>([this]() {
    return _solver().mesh().n_elements();
  }));
  _inter.variables->create<std::string>("performance_report", new Namespace::Heisenberg<std::string>([this]() {
    return _solver().stopwatch_tree().report();
  }));

  _inter.variables->create<std::string>("integrate_field", new Namespace::Heisenberg<std::string>([this]() {
    Struct_expr integrand(_vars("integrand_field").value());
    auto integral = _solver().integral_field(Qpoint_expr(integrand, _inter));
    for (unsigned i_var = 0; i_var < integrand.names.size(); ++i_var) {
      _inter.variables->assign("integral_field_" + integrand.names[i_var], integral[i_var]);
    }
    return "";
  }));
  _inter.variables->create<std::string>("integrate_surface", new Namespace::Heisenberg<std::string>([this]() {
    Struct_expr integrand(_vars("integrand_surface").value());
    auto integral = _solver().integral_surface(Boundary_expr(integrand, _inter), _solver().mesh().surface_bc_sn());
    for (unsigned i_var = 0; i_var < integrand.names.size(); ++i_var) {
      _inter.variables->assign("integral_surface_" + integrand.names[i_var], integral[i_var]);
    }
    return "";
  }));

  // load HIL code for the Case _interface
  _inter.exec(format_str(1000, "$read {%s/include/Case.hil}", config::root_dir));
  // execute input file
  _inter.exec(format_str(1000, "$read {%s}", input_file.c_str()));
}

}
