#include <filesystem>
#include <Case.hpp>
#include <Simplex_geom.hpp>
#include <read_csv.hpp>
#include <standard_atmosphere.hpp>
#include <Occt.hpp>

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
  else if (name == "nonpenetration") return new Nonpenetration;
  else if (name == "no_slip") {
    auto sub = _inter.make_sub();
    sub.exec("$thermal_bc");
    No_slip::Thermal_type therm_t;
    std::optional<double> value;
    if (sub.variables->exists("heat_flux")) { // note not recursive
      therm_t = No_slip::heat_flux;
      value = sub.variables->lookup<double>("heat_flux");
    } else if (sub.variables->exists("internal_energy")) { // note not recursive
      therm_t = No_slip::internal_energy;
      value = sub.variables->lookup<double>("internal_energy");
    } else if (sub.variables->exists("temperature")) { // note not recursive
      therm_t = No_slip::internal_energy;
      value = sub.variables->lookup<double>("temperature");
      HEXED_ASSERT(value, "thermal BC specification not understood");
      printf("temp: %e\n", value.value());
      *value *= constants::specific_gas_air/(heat_rat - 1.);
    }
    HEXED_ASSERT(value, "thermal BC specification not understood");
    return new No_slip(therm_t, value.value());
  }
  else HEXED_ASSERT(false, format_str(1000, "unrecognized boundary condition type `%s`", name.c_str()));
  return nullptr; // will never happen. just to shut up GCC warning
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
    } else if (_vars("transport_model").value() == "best") {
      visc_model.reset(new Transport_model(air_sutherland_dyn_visc));
      therm_model.reset(new Transport_model(air_sutherland_therm_cond));
    } else {
      HEXED_ASSERT(false, "unrecognized transport model specification");
    }
    _solver_ptr.reset(new Solver(*n_dim, *row_size, root_sz, _vari("local_time").value(), *visc_model, *therm_model));
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
    crit_code.push_back(_vars("surface_refine").value());
    crit_code.push_back(_vars("surface_unrefine").value());
    std::vector<std::function<bool(Element&)>> crits;
    for (std::string code : crit_code) {
      crits.emplace_back([this, code](Element& elem) {
        auto sub = _inter.make_sub();
        sub.variables->assign("is_extruded", int(!elem.tree));
        sub.variables->assign("ref_level", elem.refinement_level());
        sub.variables->assign("nom_sz", elem.nominal_size());
        sub.variables->assign("resolution_badness", elem.resolution_badness);
        auto params = elem.storage_params();
        Eigen::Vector3d center;
        center.setZero();
        for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
          center += elem.vertex(i_vert).pos;
        }
        center /= params.n_vertices();
        for (int i_dim = 0; i_dim < 3; ++i_dim) {
          sub.variables->assign("center" + std::to_string(i_dim), center(i_dim));
        }
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

  _inter.variables->create<int>("make_layers", new Namespace::Heisenberg<int>([this]() {
    _solver().mesh().disconnect_boundary(_solver().mesh().surface_bc_sn());
    _solver().mesh().extrude(Layer_sequence(_vard("wall_spacing").value(), 1));
    _solver().mesh().connect_rest(_solver().mesh().surface_bc_sn());
    _solver().calc_jacobian();
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
    std::string wd = _vars("working_dir").value();
    std::string suffix = _vars("vis_file_suffix").value();
    if (_vari("vis_field").value()) {
      State_variables sv;
      Art_visc_coef avc;
      std::vector<const Qpoint_func*> to_vis;
      to_vis.push_back(&sv);
      if (_vard("art_visc_constant").value() > 0 || _vard("art_visc_width").value() > 0) to_vis.push_back(&avc);
      _solver().visualize_field_tecplot(Qf_concat(to_vis), wd + "field" + suffix);
    }
    if (_vari("vis_surface").value() && _has_geom) {
      _solver().visualize_surface_tecplot(_solver().mesh().surface_bc_sn(), wd + "surface" + suffix);
    }
    return 0;
  }));

  _inter.variables->create<int>("update", new Namespace::Heisenberg<int>([this]() {
    if (_vard("art_visc_width").value() > 0) {
      _solver().set_art_visc_smoothness(_vard("art_visc_width").value());
    } else if (_vard("art_visc_constant").value() > 0) {
      _solver().set_art_visc_smoothness(_vard("art_visc_constant").value());
    }
    _solver().update(_vard("max_safety").value(), _vard("max_time_step").value());
    return 0;
  }));
  _inter.variables->create<int>("n_elements", new Namespace::Heisenberg<int>([this]() {
    return _solver().mesh().n_elements();
  }));
  _inter.variables->create<std::string>("header", new Namespace::Heisenberg<std::string>([this]() {
    return _solver().iteration_status().header();
  }));
  _inter.variables->create<std::string>("report", new Namespace::Heisenberg<std::string>([this]() {
    std::string rpt = _solver().iteration_status().report();
    _solver().reset_counters();
    return rpt;
  }));
  _inter.variables->create<std::string>("performance_report", new Namespace::Heisenberg<std::string>([this]() {
    return _solver().stopwatch_tree().report();
  }));

  // load HIL code for the Case _interface
  _inter.exec(format_str(1000, "$read {%s/include/Case.hil}", config::root_dir));
  // execute input file
  _inter.exec(format_str(1000, "$read {%s}", input_file.c_str()));
}

}
