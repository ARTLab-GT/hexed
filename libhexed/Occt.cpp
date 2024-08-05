#include <hexed/Occt.hpp>
#if HEXED_USE_OCCT

#include <filesystem>
// geometry
#include <TopoDS_Iterator.hxx>
#include <BRep_Tool.hxx>
#include <BRepTools.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Geom2dAPI_ProjectPointOnCurve.hxx>
#include <Geom2dAPI_InterCurveCurve.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GC_MakeLine.hxx>
#include <GC_MakePlane.hxx>
#include <GCE2d_MakeLine.hxx>
// messages
#include <Message.hxx>
#include <Message_PrinterToReport.hxx>
// file import
#include <IGESControl_Reader.hxx>
#include <STEPControl_Reader.hxx>
#include <RWStl.hxx>
// triangulations
#include <IMeshTools_Parameters.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#if HEXED_USE_TECPLOT
#include <hexed/Tecplot_file.hpp>
#endif
#include <hexed/Simplex_geom.hpp>
#include <hexed/utils.hpp>

namespace hexed {

bool Occt::message_set = false;

void Occt::set_message() {
  if (!message_set) {
    message_set = true;
    auto& printers = Message::DefaultMessenger()->ChangePrinters();
    printers.Clear();
    opencascade::handle<Message_PrinterToReport> printer;
    printers.Append(printer);
  }
}

// recursively iterates through a shape and all its sub-shapes and invokes `callback`
// on any shape of the specified type
void iterate(const TopoDS_Shape& shape, TopAbs_ShapeEnum shape_type,
             std::function<void(const TopoDS_Shape&)> callback) {
  if (shape.ShapeType() == shape_type) callback(shape);
  for (TopoDS_Iterator it(shape); it.More(); it.Next()) iterate(it.Value(), shape_type, callback);
}

void collect_curves(std::vector<opencascade::handle<Geom2d_Curve>>& curves, const TopoDS_Shape& shape) {
  opencascade::handle<Geom_Plane> plane = GC_MakePlane(0., 0., 1., 0.); // x_2 = 0 plane
  TopLoc_Location location;
  double unused [2] {}; // used by `CurveOnPlane` to return values we don't care about
  iterate(shape, TopAbs_EDGE, [&](const TopoDS_Shape& s){
    TopoDS_Edge edge = TopoDS::Edge(s);
    curves.push_back(BRep_Tool::CurveOnPlane(edge, plane, location, unused[0], unused[1]));
  });
}

Occt::Geom::Geom(const TopoDS_Shape& shape, int n_dim, double angle, double deflection, int n_segments)
: nd{n_dim}
{
  if (nd == 2) {
    collect_curves(_curves, shape);
    _simplex2.reset(new Simplex_geom<2>(segments(shape, n_segments)));
    _simplex = _simplex2.get();
  } else if (nd == 3) {
    // collect all surfaces
    iterate(shape, TopAbs_FACE, [&](const TopoDS_Shape& s){
      TopoDS_Face face = TopoDS::Face(s);
      _surfaces.push_back(BRep_Tool::Surface(face));
      if (!(_surfaces.back()->Continuity() >= GeomAbs_C1)) {
        std::cerr << "WARNING: Surface is not C1 continuous! Falling back on triangulated intersections.\n"
                  << "If you see this warning, please contact "
                  << "Micaiah: https://artlab-gt.github.io/hexed/index.html#Contact\n";
      }
    });
    Triangulation triang = _triangulate(shape, angle, deflection);
    _simplex3.reset(new Simplex_geom<3>(triang.tris, triang.params, triang.faces));
    _simplex = _simplex3.get();
  } else throw std::runtime_error("`hexed::Occt_gom` must be either 2D or 3D.");
}

void Occt::Geom::visualize(std::string format, std::string file_name) {_simplex->visualize(format, file_name);}

Nearest_point<dyn> Occt::Geom::nearest_point(Mat<> point, double max_distance, double distance_guess) {
  HEXED_ASSERT(point.size() == nd, format_str(100, "`point` must be %iD", nd));
  return _simplex->nearest_point(point, max_distance, distance_guess);
}

std::vector<double> Occt::Geom::intersections(Mat<> point0, Mat<> point1) {
  HEXED_ASSERT(point0.size() == nd, format_str(100, "`point0` must be %iD", nd));
  HEXED_ASSERT(point1.size() == nd, format_str(100, "`point1` must be %iD", nd));
  auto inters = _simplex->simplex_intersections(point0, point1);
  if (nd == 2) return inters.points;
  Mat<3> diff = point1 - point0;
  for (unsigned i_inter = 0; i_inter < inters.points.size(); ++i_inter) {
    int i_tri = inters.inds[i_inter];
    Mat<2, 3> simplex_params = _simplex3->parameters()[i_tri];
    Mat<2> params = simplex_params(all, 0);
    for (int i_coord = 0; i_coord < 2; ++i_coord) {
      params += inters.coords[i_inter](i_coord)*(simplex_params(all, 1 + i_coord) - simplex_params(all, 0));
    }
    auto& surf = _surfaces[_simplex3->faces()[i_tri]];
    if (surf->Continuity() >= GeomAbs_C1) {
      auto error_jacobian = [point0, diff, &surf](Mat<> guess){
        Mat<dyn, dyn> err_jac(3, 4);
        gp_Pnt point;
        gp_Vec vecs[2];
        surf->D1(guess(0), guess(1), point, vecs[0], vecs[1]);
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
        Mat<3> pos {point.X(), point.Y(), point.Z()};
        pos *= 1e-3;
        err_jac(all, 0) = pos - (point0 + guess(2)*diff);
        for (int i_param = 0; i_param < 2; ++i_param) {
          Mat<3> deriv {vecs[i_param].X(), vecs[i_param].Y(), vecs[i_param].Z()};
          err_jac(all, 1 + i_param) = deriv*1e-3;
        };
        err_jac(all, 3) = -diff;
        #pragma GCC diagnostic pop
        return err_jac;
      };
      Mat<> soln = math::newton(error_jacobian, Mat<3>{params[0], params[1], inters.points[i_inter]},
                                {.ftol = 1e-10*diff.norm(), .max_iters = 100});
      if (std::abs(soln(2) - inters.points[i_inter]) < 0.5) inters.points[i_inter] = soln(2);
    }
  }
  return inters.points;
}

template <typename reader_t>
TopoDS_Shape Occt::execute_reader(std::string file_name) {
  HEXED_ASSERT(std::filesystem::exists(file_name), format_str(1000, "could not find file `%s`", file_name.c_str()));
  set_message();
  reader_t reader;
  auto result = reader.ReadFile(file_name.c_str());
  HEXED_ASSERT(result == IFSelect_RetDone, format_str(1000, "failed to read geometry file (return value %i of options "
                                                            "{Void, Done, Error, Fail, Stop})", result));
  reader.TransferRoots();
  return reader.OneShape();
}

TopoDS_Shape Occt::read(std::string file_name) {
  std::string ext = file_extension(file_name);
  if      (ext == "igs" || ext == "iges") return execute_reader<IGESControl_Reader>(file_name);
  else if (ext == "stp" || ext == "step") return execute_reader<STEPControl_Reader>(file_name);
  HEXED_THROW(format_str(1000, "`hexed::Occt::read` failed to recognize file exteinsion `.%s` (case insensitive).",
                         ext.c_str()));
  throw; // just to shut up the `-Wreturn-type`---the above statement throws the actual exception
}

Occt::Triangulation Occt::_triangulate(opencascade::handle<Poly_Triangulation> poly, int face) {
  HEXED_ASSERT(!poly.IsNull(), "handle is null");
  Triangulation triang;
  bool has_params = poly->HasUVNodes();
  for (int i_tri = 0; i_tri < poly->NbTriangles(); ++i_tri) {
    auto& triangle = poly->Triangle(i_tri + 1);
    Mat<3, 3> sim;
    for (int i_vert = 0; i_vert < 3; ++i_vert) {
      gp_Pnt point = poly->Node(triangle.Value(i_vert + 1));
      sim(all, i_vert) << point.X(), point.Y(), point.Z();
    }
    triang.tris.push_back(sim*1e-3);
    if (has_params) {
      Mat<2, 3> params;
      for (int i_vert = 0; i_vert < 3; ++i_vert) {
        gp_Pnt2d point = poly->UVNode(triangle.Value(i_vert + 1));
        params(all, i_vert) << point.X(), point.Y();
      }
      triang.params.push_back(params);
      triang.faces.push_back(face);
    }
  }
  return triang;
}

Occt::Triangulation Occt::_triangulate(TopoDS_Shape shape, double angle, double deflection) {
  // setup parameters
  IMeshTools_Parameters params;
  params.Deflection               = deflection*1e3;
  params.DeflectionInterior       = deflection*1e3;
  params.Angle                    = angle;
  params.AngleInterior            = angle;
  params.Relative                 = false;
  params.InParallel               = true;
  params.MinSize                  = Precision::Confusion();
  params.InternalVerticesMode     = false;
  params.ControlSurfaceDeflection = true;
  // generate mesh
  BRepTools::Clean(shape); // get rid of any existing triangulations
  BRepMesh_IncrementalMesh mesher(shape, params);
  // fetch triangles
  Triangulation triang;
  int i_face = 0;
  iterate(shape, TopAbs_FACE, [&](const TopoDS_Shape& s){
    TopoDS_Face face = TopoDS::Face(s);
    TopLoc_Location location;
    auto face_tri = _triangulate(BRep_Tool::Triangulation(face, location), i_face++);
    triang.tris.insert(triang.tris.end(), face_tri.tris.begin(), face_tri.tris.end());
    triang.params.insert(triang.params.end(), face_tri.params.begin(), face_tri.params.end());
    triang.faces.insert(triang.faces.end(), face_tri.faces.begin(), face_tri.faces.end());
  });
  BRepTools::Clean(shape); // get rid of triangulation to avoid messing other things up
  return triang;
}

std::vector<Mat<3, 3>> Occt::triangles(opencascade::handle<Poly_Triangulation> poly) {
  return _triangulate(poly, 0).tris;
}

std::vector<Mat<3, 3>> Occt::triangles(TopoDS_Shape shape, double angle, double deflection) {
  // setup parameters
  IMeshTools_Parameters params;
  params.Deflection               = deflection*1e3;
  params.DeflectionInterior       = deflection*1e3;
  params.Angle                    = angle;
  params.AngleInterior            = angle;
  params.Relative                 = false;
  params.InParallel               = true;
  params.MinSize                  = Precision::Confusion();
  params.InternalVerticesMode     = false;
  params.ControlSurfaceDeflection = true;
  // generate mesh
  BRepTools::Clean(shape); // get rid of any existing triangulations
  BRepMesh_IncrementalMesh mesher(shape, params);
  // fetch triangles
  std::vector<Mat<3, 3>> tris;
  iterate(shape, TopAbs_FACE, [&](const TopoDS_Shape& s){
    TopoDS_Face face = TopoDS::Face(s);
    TopLoc_Location location;
    auto face_tris = triangles(BRep_Tool::Triangulation(face, location));
    tris.insert(tris.end(), face_tris.begin(), face_tris.end());
  });
  BRepTools::Clean(shape); // get rid of triangulation to avoid messing other things up
  return tris;
}

Simplex_geom<3> Occt::triangulate(TopoDS_Shape shape, double angle, double deflection, int n_div_edge) {
  Simplex_geom<3> geom(triangles(shape, angle, deflection));
  iterate(shape, TopAbs_EDGE, [&](const TopoDS_Shape& s) {
    TopoDS_Edge edge = TopoDS::Edge(s);
    double param_bounds [2];
    auto curve = BRep_Tool::Curve(edge, param_bounds[0], param_bounds[1]);
    Array<double> points({n_div_edge + 1, 3});
    for (int i_point = 0; i_point < n_div_edge + 1; ++i_point) {
      gp_Pnt point = curve->Value(param_bounds[0] + i_point*(param_bounds[1] - param_bounds[0])/n_div_edge);
      points(i_point).vector() << point.X(), point.Y(), point.Z();
    }
    points *= 1e-3; // convert to m
    geom.add_edge(points);
  });
  return geom;
}

opencascade::handle<Poly_Triangulation> Occt::read_stl(std::string file_name, double scale) {
  set_message();
  opencascade::handle<Poly_Triangulation> poly = RWStl::ReadFile(file_name.c_str());
  HEXED_ASSERT(!poly.IsNull(), "`hexed::read_stl` failed (sorry, that's all i know)");
  // convert coordinates to mm
  for (int i_node = 0; i_node < poly->NbNodes(); ++i_node) {
    poly->SetNode(i_node + 1, poly->Node(i_node + 1).XYZ()*scale*1e3);
  }
  return poly;
}

std::vector<Mat<2, 2>> Occt::segments(const TopoDS_Shape& shape, int n_segments) {
  std::vector<Mat<2, 2>> segs;
  std::vector<opencascade::handle<Geom2d_Curve>> curves;
  collect_curves(curves, shape);
  for (auto& curve : curves) {
    Mat<2> vec;
    gp_Pnt2d pnt = curve->Value(curve->FirstParameter());
    vec << pnt.X(), pnt.Y();
    for (int i_seg = 0; i_seg < n_segments; ++i_seg) {
      Mat<2, 2> seg;
      seg(all, 0) = vec;
      double interp = (i_seg + 1.)/n_segments;
      pnt = curve->Value((1 - interp)*curve->FirstParameter() + interp*curve->LastParameter());
      vec << pnt.X(), pnt.Y();
      seg(all, 1) = vec;
      seg *= 1e-3; // convert to m
      segs.push_back(seg);
    }
  }
  return segs;
}

}
#endif
