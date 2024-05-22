#include <filesystem>
#include <cctype>
#include <Case.hpp>
#include <Simplex_geom.hpp>
#include <read_csv.hpp>
#include <standard_atmosphere.hpp>
#include <Occt.hpp>
#include <hil_properties.hpp>
#include <Csv.hpp>

namespace hexed
{

const double heat_rat = 1.4;

Solver& Case::_solver()
{
  HEXED_ASSERT(_solver_ptr, "`Solver` object does not exist", assert::User_error);
  return *_solver_ptr;
}

std::string strip_trailing_digits(std::string s)
{
  while (std::isdigit(s.back())) s.pop_back();
  return s;
}

        int Case::_vari(std::string name) {return _inter.variables->get<        int>(name);}
     double Case::_vard(std::string name) {return _inter.variables->get<     double>(name);}
std::string Case::_vars(std::string name) {return _inter.variables->get<std::string>(name);}

Mat<> Case::_get_vector(std::string name, int size)
{
  Mat<> vec(size);
  for (int i = 0; i < size; ++i) {
    HEXED_ASSERT(_inter.variables->lookup<double>(name + std::to_string(i)), "must specify all components of `" + name + "` or none", assert::User_error);
    vec(i) = _vard(name + std::to_string(i));
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
  Mat<> freestream = _get_vector("freestream", _vari("n_dim") + 2);
  if      (name == "characteristic") return new Riemann_invariants(freestream);
  else if (name == "freestream") return new Freestream(freestream);
  else if (name == "pressure_outflow") return new Pressure_outflow(_vard("freestream_pressure"));
  else if (name == "outflow") return new Outflow;
  else if (name == "nonpenetration") return new Nonpenetration;
  else if (name == "no_slip") {
    auto sub = _inter.make_sub();
    sub.exec("$thermal_bc");
    std::shared_ptr<Thermal_bc> thermal;
    if (sub.variables->exists("heat_flux")) { // note not recursive
      thermal = std::make_shared<Prescribed_heat_flux>(sub.variables->lookup<double>("heat_flux").value());
    } else if (sub.variables->exists("emissivity") || sub.variables->exists("heat_transfer_coef")) {
      auto equilibrium = std::make_shared<Thermal_equilibrium>();
      HEXED_ASSERT(sub.variables->exists("heat_transfer_coef") == sub.variables->exists("temperature"),
        "must specify both surface heat_transfer_coef and temperature or neither", assert::User_error);
      if (sub.variables->exists("emissivity")) equilibrium->emissivity = sub.variables->lookup<double>("emissivity").value();
      if (sub.variables->exists("heat_transfer_coef")) {
        equilibrium->heat_transfer_coef = sub.variables->lookup<double>("heat_transfer_coef").value();
        equilibrium->temperature = sub.variables->lookup<double>("temperature").value();
      }
      thermal = equilibrium;
    } else if (sub.variables->exists("internal_energy")) { // note not recursive
      thermal = std::make_shared<Prescribed_energy>(sub.variables->lookup<double>("internal_energy").value());
    } else if (sub.variables->exists("temperature")) { // note not recursive
      double energy = sub.variables->lookup<double>("temperature").value()*constants::specific_gas_air/(heat_rat - 1.);
      thermal = std::make_shared<Prescribed_energy>(energy);
    }
    HEXED_ASSERT(thermal, "thermal BC specification not understood", assert::User_error);
    return new No_slip(thermal, _vard("heat_flux_coercion"));
  }
  else HEXED_ASSERT(false, format_str(1000, "unrecognized boundary condition type `%s`", name.c_str()), assert::User_error);
  return nullptr; // will never happen. just to shut up GCC warning
}

std::string Case::_iteration_suffix()
{
  return format_str(100, "iter%.*i", _vari("iter_width"), _vari("iteration"));
}

void force_symlink(const std::filesystem::path& target, const std::filesystem::path& link)
{
  if (std::filesystem::exists(link)) std::filesystem::remove(link);
  std::filesystem::create_symlink(target, link);
}

std::vector<Flow_bc*> Case::_make_extremal_bcs()
{
  std::vector<Flow_bc*> bcs;
  for (int i_dim = 0; i_dim < _vari("n_dim"); ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      bcs.push_back(_make_bc(_vars(format_str(50, "extremal_bc%i%i", i_dim, sign))));
    }
  }
  return bcs;
}

Surface_geom* Case::_make_geom()
{
  int nd = _vari("n_dim");
  std::vector<Surface_geom*> geoms;
  for (int i_geom = 0;; ++i_geom) {
    auto geom = _inter.variables->lookup<std::string>("geom" + std::to_string(i_geom));
    if (!geom) break;
    unsigned dot = geom->rfind('.');
    HEXED_ASSERT(dot < geom->size(), "file name must contain extension to infer format", assert::User_error);
    HEXED_ASSERT(std::filesystem::exists(geom.value()), format_str(1000, "geometry file `%s` not found", geom->c_str()), assert::User_error);
    std::string case_sensitive(geom->begin() + dot + 1, geom->end());
    std::string ext = case_sensitive;
    for (char& c : ext) c = tolower(c);
    if (ext == "csv") {
      HEXED_ASSERT(nd == 2, "3D geometry in CSV format is not supported", assert::User_error);
      auto data = read_csv(*geom);
      HEXED_ASSERT(data.cols() >= nd, "CSV geometry file must have at least n_dim columns", assert::User_error);
      geoms.emplace_back(new Simplex_geom<2>(segments(data.transpose())));
    #if HEXED_USE_OCCT
    } else if (ext == "igs" || ext == "iges" || ext == "stp" || ext == "step") {
      auto shape = Occt::read(*geom);
      if (nd == 2) {
        geoms.emplace_back(new Simplex_geom<2>(Occt::segments(shape, _vari("geom_n_segments"))));
      } else if (nd == 3) {
        auto ptr = new Simplex_geom<3>(Occt::triangles(shape, _vard("max_angle"), _vard("max_deflection")));
        std::string vis_name = format_str(1000, "%sgeom%i_triangulation", _vars("working_dir").c_str(), i_geom);
        #if HEXED_USE_XDMF
        ptr->visualize("xdmf", vis_name);
        #elif HEXED_USE_TECPLOT
        ptr->visualize("tecplot", vis_name);
        #endif
        geoms.emplace_back(ptr);
      }
    } else if (ext == "stl") {
      HEXED_ASSERT(nd == 3, "STL format is only supported for 3D", assert::User_error);
      geoms.emplace_back(new Simplex_geom<3>(Occt::triangles(Occt::read_stl(geom.value()))));
    #endif
    } else {
      HEXED_ASSERT(false, format_str(1000, "file extension `%s` not recognized", case_sensitive.c_str()), assert::User_error);
    }
  }
  return geoms.empty() ? nullptr : new Compound_geom(geoms);
  return nullptr;
}

std::string Case::_assignment(std::string var_name)
{
  std::string statement = var_name + " = ";
  if      (_inter.variables->lookup<        int>(var_name)) statement += std::to_string(_vari(var_name));
  else if (_inter.variables->lookup<     double>(var_name)) statement += format_str(100, "%.20e", _vard(var_name));
  else if (_inter.variables->lookup<std::string>(var_name)) statement += "{" + _vars(var_name) + "}";
  return statement;
}

Case::Case(std::string input_script)
: _printers{std::make_shared<Printer_set>()}
{
  _inter.variables->assign("input_script", input_script);
  _inter.variables->assign("version_major", config::version_major);
  _inter.variables->assign("version_minor", config::version_minor);
  _inter.variables->assign("version_patch", config::version_patch);
  _inter.variables->assign<std::string>("commit", config::commit);

  // create custom Heisenberg variables

  _inter.variables->create("setup_output", new Namespace::Heisenberg<std::string>([this]() {
    _output_file.reset(new std::ofstream(_vars("working_dir") + "output.txt"));
    for (auto* printer : {&_printers->info, &_printers->warn, &_printers->error}) {
      printer->printers.emplace_back(std::make_shared<Stream_printer>(*_output_file));
    }
    _inter.printer = _printers;
    char utc [100];
    std::time_t time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::strftime(utc, 100, "%Y-%m-%d %H:%M:%S", std::gmtime(&time));
    _printers->info(format_str(1000, "Commencing simulation with Hexed version %i.%i.%i (commit %s) at %s UTC (%i Unix Time).\n",
                               config::version_major, config::version_minor, config::version_patch, config::commit.c_str(), utc, time));
    return "";
  }));

  _inter.variables->create("setup_parameters", new Namespace::Heisenberg<std::string>([this]() {
    // setup storage parameters
    auto n_dim = _inter.variables->lookup<int>("n_dim");
    HEXED_ASSERT(n_dim && (n_dim.value() > 0) && (n_dim.value() <= 3),
                 "`n_dim` must be defined as an integer in [1, 3]", assert::User_error);
    auto row_size = _inter.variables->lookup<int>("row_size");
    HEXED_ASSERT(row_size.value() >= 2 && row_size.value() <= config::max_row_size,
                 format_str(300, "`row_size` must be between 2 and %i", config::max_row_size), assert::User_error);
    // compute freestream
    Mat<> freestream(*n_dim + 2);
    if (_inter.variables->lookup<double>("freestream0")) freestream = _get_vector("freestream", *n_dim + 2);
    else {
      if (_inter.variables->lookup<double>("altitude")) {
        HEXED_ASSERT(!_inter.variables->lookup<double>("freestream_temperature"), "cannot specify both altitude and temperature (consider `temperature_offset`)",
                      assert::User_error);
        auto dens_pres = standard_atmosphere(_vard("altitude"), _vard("temperature_offset"));
        _inter.variables->assign<double>("freestream_density", dens_pres[0]);
        _inter.variables->assign<double>("freestream_pressure", dens_pres[1]);
      }
      HEXED_ASSERT(  _inter.variables->lookup<double>("freestream_density").has_value()
                   + _inter.variables->lookup<double>("freestream_pressure").has_value()
                   + _inter.variables->lookup<double>("freestream_temperature").has_value() == 2,
                   "exactly two of freestream density, pressure, and temperature must be specified", assert::User_error);
      if (_inter.variables->lookup<double>("freestream_density")) {
        freestream(*n_dim) = _vard("freestream_density");
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
      auto direction = full_direction(Eigen::seqN(0, *n_dim));
      if (_inter.variables->lookup<double>("freestream_velocity0")) {
        veloc = _get_vector("freestream_velocity", *n_dim);
        direction = veloc.normalized();
      } else {
        _inter.variables->assign<double>("freestream_sound_speed", std::sqrt(heat_rat*constants::specific_gas_air*_vard("freestream_temperature")));
        if (_inter.variables->lookup<double>("freestream_speed")) _inter.variables->assign<double>("freestream_mach", _vard("freestream_speed")/_vard("freestream_sound_speed"));
        else _inter.variables->assign<double>("freestream_speed", _vard("freestream_mach")*_vard("freestream_sound_speed"));
        if (_inter.variables->lookup<double>("freestream_direction0")) direction = _get_vector("freestream_direction", *n_dim).normalized();
        else {
          direction.setUnit(*n_dim, 0);
          if (*n_dim == 2) {
            direction = Eigen::Rotation2D<double>(_vard("attack"))*direction;
          }
          if (*n_dim == 3) {
            direction = Eigen::AngleAxis<double>(-_vard("attack"  ), Eigen::Vector3d::Unit(1))*direction;
            direction = Eigen::AngleAxis<double>( _vard("sideslip"), Eigen::Vector3d::Unit(2))*direction;
          }
        }
        veloc = _vard("freestream_speed")*direction;
        _set_vector("freestream_velocity", veloc);
      }
      _set_vector("freestream_direction", full_direction);
      freestream(Eigen::seqN(0, *n_dim)) = _vard("freestream_density")*veloc;
      freestream(*n_dim) = _vard("freestream_density");
      freestream(*n_dim + 1) = _vard("freestream_pressure")/(heat_rat - 1) + .5*_vard("freestream_density")*veloc.squaredNorm();
      freestream.conservativeResize(5);
      freestream(Eigen::seqN(*n_dim + 2, 5 - (*n_dim + 2))).setZero();
      _set_vector("freestream", freestream);
    }
    // create history monitors
    _monitor_expr.reset(new Struct_expr(_vars("monitor_vars")));
    for (std::string name : _monitor_expr->names) {
      _monitors.emplace_back(_vard("monitor_window"), _vari("monitor_samples"));
      _inter.variables->assign(name + "_min", -huge);
      _inter.variables->assign(name + "_max",  huge);
    }
    return "";
  }));

  _inter.variables->create("create_solver", new Namespace::Heisenberg<std::string>([this]() {
    int n_dim = _vari("n_dim");
    // evaluate dimensions
    Mat<dyn, dyn> mesh_extremes(n_dim, 2);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      for (int sign = 0; sign < 2; ++sign) {
        mesh_extremes(i_dim, sign) = _vard(format_str(50, "mesh_extreme%i%i", i_dim, sign));
      }
    }
    HEXED_ASSERT((mesh_extremes(all, 1) - mesh_extremes(all, 0)).minCoeff() > 0, "all mesh dimensions must be positive!", assert::User_error);
    double root_size = (mesh_extremes(all, 1) - mesh_extremes(all, 0)).maxCoeff();
    // construct molecular transport models
    std::vector<std::string> transport_phenomena {"viscosity", "conductivity"};
    std::vector<Transport_model> transport_models;
    for (std::string name : transport_phenomena) {
      auto sub = _inter.make_sub();
      Struct_expr model(_vars(name + "_model"));
      model.eval(sub);
      if (model.names.empty()) transport_models.emplace_back(inviscid);
      else if (sub.variables->exists("offset")) {
        transport_models.emplace_back(Transport_model::sutherland(sub.variables->lookup<double>("ref_value").value(),
                                                                  sub.variables->lookup<double>("ref_temperature").value(),
                                                                  sub.variables->lookup<double>("offset").value()));
      } else HEXED_ASSERT(false, format_str(200, "invalid transport model specification for %s", name), assert::User_error);
    }
    // setup actual solver
    _solver_ptr.reset(new Solver(n_dim, _vari("row_size"), root_size, true, transport_models[0], transport_models[1], _inter.variables, _printers));
    _solver().mesh().add_tree(_make_extremal_bcs(), mesh_extremes(all, 0));
    _solver().set_fix_admissibility(_vari("fix_therm_admis"));
    return "";
  }));

  _inter.variables->create("init_refinement", new Namespace::Heisenberg<std::string>([this]() {
    for (int i = 0; i < _vari("init_ref_level"); ++i) _solver().mesh().update();
    _solver().calc_jacobian();
    return "";
  }));

  _inter.variables->create("add_geom", new Namespace::Heisenberg<std::string>([this]() {
    Surface_geom* geom = _make_geom();
    if (geom) {
      _has_geom = true;
      _solver().mesh().set_surface(geom, _make_bc(_vars("surface_bc")), _get_vector("flood_fill_start", _vari("n_dim")));
      for (int i_smooth = 0; i_smooth < _vari("n_smooth"); ++i_smooth) _solver().mesh().relax(0.5);
      _solver().calc_jacobian();
    }
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
        hil_properties::element(*sub.variables, elem);
        sub.exec(code);
        return sub.variables->lookup<int>("return").value();
      });
    }
    Jac_inv_det_func jidf;
    _solver().set_uncertainty(Elem_nonsmooth(jidf));
    _solver().mesh().set_unref_locks(criteria::if_extruded);
    bool changed = _solver().mesh().update(crits[0], crits[1]);
    for (int i_smooth = 0; i_smooth < _vari("n_smooth"); ++i_smooth) _solver().mesh().relax(0.5);
    _solver().calc_jacobian();
    return changed;
  }));

  _inter.variables->create("split_layers", new Namespace::Heisenberg<std::string>([this]() {
    _solver().mesh().disconnect_boundary(_solver().mesh().surface_bc_sn());
    auto sub = _inter.make_sub();
    std::vector<double> split_points = Struct_expr(_vars("layer_split_points")).eval(sub);
    double prev_split = 1.;
    for (double split : split_points) {
      _solver().mesh().extrude(true, split/prev_split, true);
      prev_split = split;
    }
    _solver().mesh().connect_rest(_solver().mesh().surface_bc_sn());
    _solver().calc_jacobian();
    return "";
  }));

  _inter.variables->create("init_state", new Namespace::Heisenberg<std::string>([this]() {
    #if HEXED_OBSESSIVE_TIMING
    int nd = _inter.variables->lookup<int>("n_dim").value();
    if (nd == 2) _inter.printer->print(Simplex_geom<2>::performance_report());
    if (nd == 3) _inter.printer->print(Simplex_geom<3>::performance_report());
    #endif
    _solver().initialize(Spacetime_expr(Struct_expr(_vars("init_cond")), _inter));
    return "";
  }));

  _inter.variables->create("read_mesh", new Namespace::Heisenberg<std::string>([this]() {
    _printers->info("reading mesh... ");
    Surface_geom* geom = _make_geom();
    _solver().read_mesh(_vars("input_data"), _make_extremal_bcs(), geom, geom ? _make_bc(_vars("surface_bc")) : nullptr);
    _printers->info("done\n");
    return "";
  }));
  _inter.variables->create("read_state", new Namespace::Heisenberg<std::string>([this]() {
    _printers->info("reading state... ");
    _solver().read_state(_vars("input_data"));
    _printers->info("done\n");
    return "";
  }));
  _inter.variables->create("write_mesh", new Namespace::Heisenberg<std::string>([this]() {
    _printers->info("writing mesh... ");
    std::string file_name = _vars("working_dir") + _iteration_suffix();
    _solver().mesh().write(file_name);
    force_symlink(_iteration_suffix() + ".mesh.h5", _vars("working_dir") + "latest.mesh.h5");
    _printers->info("done\n");
    return "";
  }));
  _inter.variables->create("write_state", new Namespace::Heisenberg<std::string>([this]() {
    _printers->info("writing state... ");
    std::string file_name = _vars("working_dir") + _iteration_suffix();
    _solver().write_state(file_name);
    force_symlink(_iteration_suffix() + ".state.h5", _vars("working_dir") + "latest.state.h5");
    _printers->info("done\n");
    return "";
  }));
  _inter.variables->create("write_status", new Namespace::Heisenberg<std::string>([this]() {
    _printers->info("writing status... ");
    std::ofstream status_file(_vars("working_dir") + _iteration_suffix() + ".status.hil");
    std::vector<std::string> no_write {"working_dir", "input_data"};
    for (std::string name : _inter.variables->names()) {
      if (name.substr(0, 6) != "hexed_" && std::none_of(no_write.begin(), no_write.end(), [name](std::string nw){return name == nw;})) {
        status_file << _assignment(name) + "\n";
      }
    }
    status_file.close();
    force_symlink(_iteration_suffix() + ".status.hil", _vars("working_dir") + "latest.status.hil");
    _printers->info("done\n");
    return "";
  }));

  _inter.variables->create("visualize", new Namespace::Heisenberg<std::string>([this]() {
    _printers->info("visualizing... ");
    std::string wd = _vars("working_dir");
    std::string suffix = "_" + _iteration_suffix();
    int n_sample = _vari("vis_n_sample");
    std::vector<std::string> vis_objects {"surface", "field"};
    for (int i_contour = 0; ; ++i_contour) {
      std::string name = "contour" + std::to_string(i_contour);
      if (_inter.variables->lookup<std::string>(name)) vis_objects.push_back(name);
      else break;
    }
    for (std::string v : vis_objects) if (_vari("vis_" + strip_trailing_digits(v))) {
      for (std::string format : {"xdmf", "tecplot", "csv"}) if (_vari("vis_" + format)) {
        Struct_expr vis_vars(_vars("vis_" + strip_trailing_digits(v) + "_vars"));
        for (bool edges : {false, true}) {
          std::string name = v;
          if (edges) name = name + "_edges";
          if (!edges || (_vari("vis_edges") && _vari("n_dim") + (v == "field") > 2)) {
            std::string file_name = wd + name + suffix;
            if (v == "surface") {
              _solver().visualize_surface(format, file_name, _solver().mesh().surface_bc_sn(), Boundary_expr(vis_vars, _inter), n_sample, edges);
            } else if (v == "field") {
              _solver().visualize_field(format, file_name, Qpoint_expr(vis_vars, _inter), n_sample, edges);
              if (_vari("vis_skew")) _solver().visualize_field(format, wd + "skew" + suffix, Equiangle_skewness(), n_sample, edges);
              if (_vari("vis_lts_constraints")) _solver().vis_lts_constraints(format, wd + "lts_constraints" + suffix, n_sample);
            } else if (!edges) { // vis_type == contour0, contour1, etc
              std::string contour_expr = _vars("vis_contour_vars") + v + "_var = " + _vars(v) + ";";
              _solver().visualize_contour(format, file_name, Qpoint_expr(contour_expr, _inter), Qpoint_expr(vis_vars, _inter), n_sample);
            }
            if (format == "xdmf") {
              std::string latest = wd + name + "_latest1.xmf";
              if (std::filesystem::exists(latest)) {
                std::filesystem::copy_file(latest, wd + name + "_latest0.xmf", std::filesystem::copy_options::overwrite_existing);
              }
              if (std::filesystem::exists(file_name + ".xmf")) {
                std::filesystem::copy_file(file_name + ".xmf", latest, std::filesystem::copy_options::overwrite_existing);
              }
            }
          }
        }
      }
    }
    _printers->info("done\n");
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
    Struct_expr vars(_vars("print_vars"));
    for (std::string name : vars.names) {
      int width = std::max<int>(15, name.size());
      header += format_str(1000, "%*s, ", width, name.c_str());
    }
    header.erase(header.end() - 2, header.end());
    return header;
  }));

  _inter.variables->create<std::string>("compute_residuals", new Namespace::Heisenberg<std::string>([this]() {
    int nd = _solver().storage_params().n_dim;
    Physical_residual phys_resid;
    _solver().compute_residual();
    auto res = _solver().integral_field(Pow(phys_resid, 2));
    for (int i_dim = 1; i_dim < nd; ++i_dim) res[0] += res[i_dim];
    for (double& r : res) {
      r = std::sqrt(r);
      HEXED_ASSERT(!std::isnan(r), "residual is NaN", assert::Numerical_exception);
    }
    _inter.variables->assign("residual_momentum", res[0]);
    _inter.variables->assign("residual_density", res[nd]);
    _inter.variables->assign("residual_energy", res[nd + 1]);
    return "";
  }));

  _inter.variables->create<std::string>("compute_lts_constraints", new Namespace::Heisenberg<std::string>([this]() {
    _solver().compute_lts_constraints();
    return "";
  }));

  _inter.variables->create<std::string>("report", new Namespace::Heisenberg<std::string>([this]() {
    std::string report = "";
    Struct_expr vars(_vars("print_vars"));
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
    HEXED_ASSERT(_vari("mesh_init"), "attempt to update flow when mesh has not been created", assert::User_error);
    bool avw = _vard("art_visc_width") > 0;
    bool avc = _vard("art_visc_constant") > 0;
    int iter = _vari("iteration");
    int print_freq = _vari("print_freq");
    int n = iter ? print_freq - iter%print_freq : 1;
    for (int i = 0; i < n; ++i) {
      ++iter;
      if (_inter.variables->get<int>("diffusive_admissibility")) _solver().set_art_visc_admis();
      if (_inter.variables->get<int>("elementwise_art_visc")) {
        _solver().update_art_visc_elwise(_vard("art_visc_width"), _vari("elementwise_art_visc_pde"));
      } else if (avw) {
        _solver().update_art_visc_smoothness(_vard("art_visc_width"));
      } else if (avc) {
        _solver().set_art_visc_constant(_vard("art_visc_constant"));
      }
      _solver().update();
    }
    _inter.variables->assign("iteration", iter);
    auto sub = _inter.make_sub();
    auto vals = _monitor_expr->eval(sub);
    for (unsigned i_monitor = 0; i_monitor < _monitor_expr->names.size(); ++i_monitor) {
      _monitors[i_monitor].add_sample(iter, vals[i_monitor]);
      _inter.variables->assign(_monitor_expr->names[i_monitor] + "_min", _monitors[i_monitor].min());
      _inter.variables->assign(_monitor_expr->names[i_monitor] + "_max", _monitors[i_monitor].max());
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
    Struct_expr integrand(_vars("integrand_field"));
    auto integral = _solver().integral_field(Qpoint_expr(integrand, _inter));
    for (unsigned i_var = 0; i_var < integrand.names.size(); ++i_var) {
      _inter.variables->assign("integral_field_" + integrand.names[i_var], integral[i_var]);
    }
    return "";
  }));
  _inter.variables->create<std::string>("integrate_surface", new Namespace::Heisenberg<std::string>([this]() {
    Struct_expr integrand(_vars("integrand_surface"));
    auto integral = _solver().integral_surface(Boundary_expr(integrand, _inter), _solver().mesh().surface_bc_sn());
    for (unsigned i_var = 0; i_var < integrand.names.size(); ++i_var) {
      _inter.variables->assign("integral_surface_" + integrand.names[i_var], integral[i_var]);
    }
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
    _inter.exec(format_str(1000, "$read {%s}", input_script.c_str()));
  } catch (const assert::User_error& except) {
    if (_printers) _printers->error("User error: ", true);
    throw except;
  } catch (const assert::Numerical_exception& except) {
    _printers->error("Numerical exception: ", true);
    _printers->error(except.what());
    _printers->error("\nTerminating simulation.\n", true);
    if (_solver_ptr) _inter.exec("write_mesh; write_state; write_status; visualize;");
    _inter.variables->assign("failed", 1);
  }
}

}
