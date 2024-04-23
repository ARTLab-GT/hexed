#include <Occt.hpp>
#if HEXED_USE_OCCT

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
// rendering
#include <AIS_InteractiveContext.hxx>
#include <AIS_Shape.hxx>
#include <Aspect_DisplayConnection.hxx>
#include <OpenGl_GraphicDriver.hxx>
#include <V3d_View.hxx>
#include <Xw_Window.hxx>
// file import
#include <IGESControl_Reader.hxx>
#include <STEPControl_Reader.hxx>
#include <RWStl.hxx>
// triangulations
#include <IMeshTools_Parameters.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#if HEXED_USE_TECPLOT
#include <Tecplot_file.hpp>
#endif
#include <Simplex_geom.hpp>

namespace hexed
{

bool Occt::message_set = false;

void Occt::set_message()
{
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
void iterate(const TopoDS_Shape& shape, TopAbs_ShapeEnum shape_type, std::function<void(const TopoDS_Shape&)> callback)
{
  if (shape.ShapeType() == shape_type) callback(shape);
  for (TopoDS_Iterator it(shape); it.More(); it.Next()) iterate(it.Value(), shape_type, callback);
}

void collect_curves(std::vector<opencascade::handle<Geom2d_Curve>>& curves, const TopoDS_Shape& shape)
{
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
    _simplex.reset(new Simplex_geom<2>(segments(shape, n_segments)));
  } else if (nd == 3) {
    // collect all surfaces
    iterate(shape, TopAbs_FACE, [&](const TopoDS_Shape& s){
      TopoDS_Face face = TopoDS::Face(s);
      _surfaces.push_back(BRep_Tool::Surface(face));
    });
    _simplex.reset(new Simplex_geom<3>(triangles(shape, angle, deflection)));
  } else throw std::runtime_error("`hexed::Occt_gom` must be either 2D or 3D.");
}

void Occt::Geom::visualize(std::string format, std::string file_name) {_simplex->visualize(format, file_name);}

Nearest_point<dyn> Occt::Geom::nearest_point(Mat<> point, double max_distance, double distance_guess)
{
  HEXED_ASSERT(point.size() == nd, format_str(100, "`point` must be %iD", nd));
  return _simplex->nearest_point(point, max_distance, distance_guess);
}

std::vector<double> Occt::Geom::intersections(Mat<> point0, Mat<> point1)
{
  HEXED_ASSERT(point0.size() == nd, format_str(100, "`point0` must be %iD", nd));
  HEXED_ASSERT(point1.size() == nd, format_str(100, "`point1` must be %iD", nd));
  // convert to mm
  Mat<> scaled0 = 1000*point0;
  Mat<> scaled1 = 1000*point1;
  double dist = (scaled1 - scaled0).norm();
  std::vector<double> sects;
  if (nd == 2) {
    // compute the line through the given points
    gp_Pnt2d pnt0(scaled0(0), scaled0(1));
    gp_Pnt2d pnt1(scaled1(0), scaled1(1));
    opencascade::handle<Geom2d_Line> line = GCE2d_MakeLine(pnt0, pnt1);
    // iterate through _curves and compute the indersections with each
    for (auto& curve : _curves) {
      // find intersections
      Geom2dAPI_InterCurveCurve inter(line, curve);
      int n = inter.NbPoints();
      for (int i = 0; i < n; ++i) {
        // 2d intersector doesn't seem to have a `Parameters` member,
        // so we have to compute the parametric representation ourselves
        // as the least squares solution to `t*(scaled1 - scaled0) = point - scaled0`
        gp_Pnt2d occt_point = inter.Point(i + 1);
        Mat<2> point{occt_point.X(), occt_point.Y()};
        Mat<2> lhs = scaled1 - scaled0;
        Mat<2> rhs = point - scaled0;
        sects.push_back(lhs.dot(rhs)/lhs.squaredNorm());
      }
    }
  } else {
    // compute the line through the given points
    gp_Pnt pnt0(scaled0(0), scaled0(1), scaled0(2));
    gp_Pnt pnt1(scaled1(0), scaled1(1), scaled1(2));
    opencascade::handle<Geom_Line> line = GC_MakeLine(pnt0, pnt1);
    opencascade::handle<Geom_Curve> curve = line;
    // iterate through _surfaces and compute the indersections with each
    for (auto& surface : _surfaces) {
      // compute intersections
      GeomAPI_IntCS inter(curve, surface);
      HEXED_ASSERT(inter.IsDone(), "line/surface intersection failed in OCCT kernel", assert::Numerical_exception);
      // translate to our parametric format
      int n = inter.NbPoints();
      for (int i = 0; i < n; ++i) {
        double params [3];
        inter.Parameters(i + 1, params[0], params[1], params[2]);
        sects.push_back(params[2]/dist);
      }
    }
  }
  return math::correct_values(_simplex->intersections(point0, point1), sects);
}

void Occt::write_image(const TopoDS_Shape& shape, std::string file_name, Mat<3> eye_pos, Mat<3> look_at_pos, int resolution)
{
  // general setup
  opencascade::handle<Aspect_DisplayConnection> displayConnection = new Aspect_DisplayConnection();
  opencascade::handle<OpenGl_GraphicDriver> graphicDriver = new OpenGl_GraphicDriver(displayConnection);
  opencascade::handle<V3d_Viewer> viewer = new V3d_Viewer(graphicDriver);
  viewer->SetDefaultLights();
  viewer->SetLightOn();
  opencascade::handle<AIS_InteractiveContext> context = new AIS_InteractiveContext(viewer);
  opencascade::handle<V3d_View> view = viewer->CreateView();
  opencascade::handle<Xw_Window> win = new Xw_Window(graphicDriver->GetDisplayConnection(), "", 0, 0, resolution, resolution);
  win->SetVirtual(true);
  view->SetWindow(win);
  view->SetBackgroundColor(Quantity_Color(Quantity_NOC_BLACK));
  view->MustBeResized();
  view->AutoZFit();
  // add shapes
  opencascade::handle<AIS_Shape> shaded = new AIS_Shape(shape);
  context->Display(shaded, false);
  context->SetDisplayMode(shaded, AIS_Shaded, false);
  opencascade::handle<AIS_Shape> wireframe = new AIS_Shape(shape);
  context->Display(wireframe, false);
  context->SetDisplayMode(wireframe, AIS_WireFrame, false);
  view->SetFront();
  // set view orientation
  eye_pos *= 1e3;
  look_at_pos *= 1e3;
  view->SetEye(eye_pos(0), eye_pos(1), eye_pos(2));
  view->SetAt(look_at_pos(0), look_at_pos(1), look_at_pos(2));
  view->FitAll(.2);
  // add coordinate axes
  view->TriedronDisplay(Aspect_TOTP_LEFT_LOWER, Quantity_NOC_WHITE, .1);
  // render/save
  view->Redraw();
  view->Dump(file_name.c_str());
}

template <typename reader_t>
TopoDS_Shape Occt::execute_reader(std::string file_name)
{
  set_message();
  reader_t reader;
  auto result = reader.ReadFile(file_name.c_str());
  HEXED_ASSERT(result == IFSelect_RetDone, "failed to read geometry file");
  reader.TransferRoots();
  return reader.OneShape();
}

TopoDS_Shape Occt::read(std::string file_name)
{
  unsigned extension_start = file_name.find_last_of(".");
  HEXED_ASSERT(extension_start != std::string::npos, "`file_name` has no extension");
  std::string case_sensitive = file_name.substr(extension_start + 1, std::string::npos);
  std::string ext = case_sensitive;
  for (char& c : ext) c = tolower(c);
  if      (ext == "igs" || ext == "iges") return execute_reader<IGESControl_Reader>(file_name);
  else if (ext == "stp" || ext == "step") return execute_reader<STEPControl_Reader>(file_name);
  throw std::runtime_error(format_str(1000, "`hexed::Occt::read` failed to recognize file exteinsion `.%s`.", case_sensitive.c_str()));
}

std::vector<Mat<3, 3>> Occt::triangles(opencascade::handle<Poly_Triangulation> poly)
{
  HEXED_ASSERT(!poly.IsNull(), "handle is null");
  std::vector<Mat<3, 3>> sims;
  for (int i_tri = 0; i_tri < poly->NbTriangles(); ++i_tri) {
    auto& triangle = poly->Triangle(i_tri + 1);
    Mat<3, 3> sim;
    for (int i_vert = 0; i_vert < 3; ++i_vert) {
      gp_Pnt point = poly->Node(triangle.Value(i_vert + 1));
      sim(all, i_vert) << point.X(), point.Y(), point.Z();
    }
    sims.push_back(sim*1e-3);
  }
  return sims;
}

std::vector<Mat<3, 3>> Occt::triangles(TopoDS_Shape shape, double angle, double deflection)
{
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

opencascade::handle<Poly_Triangulation> Occt::read_stl(std::string file_name, double scale)
{
  set_message();
  opencascade::handle<Poly_Triangulation> poly = RWStl::ReadFile(file_name.c_str());
  HEXED_ASSERT(!poly.IsNull(), "`hexed::read_stl` failed (sorry, that's all i know)");
  // convert coordinates to mm
  for (int i_node = 0; i_node < poly->NbNodes(); ++i_node) {
    poly->SetNode(i_node + 1, poly->Node(i_node + 1).XYZ()*scale*1e3);
  }
  return poly;
}

std::vector<Mat<2, 2>> Occt::segments(const TopoDS_Shape& shape, int n_segments)
{
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
