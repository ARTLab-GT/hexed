#include <Occt_geom.hpp>
#if HEXED_USE_OCCT

// geometry
#include <TopoDS_Iterator.hxx>
#include <BRep_Tool.hxx>
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
#if HEXED_USE_TECPLOT
#include <Tecplot_file.hpp>
#endif
#include <iostream>

namespace hexed
{

bool Occt_geom::message_set = false;

void Occt_geom::set_message()
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

Occt_geom::Occt_geom(const TopoDS_Shape& shape, int n_dim)
: nd{n_dim}
{
  if (nd == 2) {
    collect_curves(curves, shape);
  } else if (nd == 3) {
    // collect all surfaces
    iterate(shape, TopAbs_FACE, [&](const TopoDS_Shape& s){
      TopoDS_Face face = TopoDS::Face(s);
      surfaces.push_back(BRep_Tool::Surface(face));
    });
  } else throw std::runtime_error("`hexed::Occt_gom` must be either 2D or 3D.");
}

void Occt_geom::visualize(std::string file_name)
{
  int n_div = 20;
  if (nd == 3) {
    Tecplot_file file(file_name, 3, {"real"}, 0.);
    for (unsigned i_surf = 0; i_surf < surfaces.size(); ++i_surf) {
      auto& surf = surfaces[i_surf];
      double param_bounds [2][2];
      surf->Bounds(param_bounds[0][0], param_bounds[0][1], param_bounds[1][0], param_bounds[1][1]);
      Mat<dyn, dyn> data(math::pow(n_div + 1, 2), 3);
      Mat<dyn> real(math::pow(n_div + 1, 2));
      int node_coords [2];
      for (node_coords[0] = 0; node_coords[0] <= n_div; ++node_coords[0]) {
        for (node_coords[1] = 0; node_coords[1] <= n_div; ++node_coords[1]) {
          double params [2];
          for (int i_dim = 0; i_dim < 2; ++i_dim) {
            double interp = double(node_coords[i_dim])/n_div;
            params[i_dim] = (1 - interp)*param_bounds[i_dim][0] + interp*param_bounds[i_dim][1];
          }
          auto pnt = surf->Value(params[0], params[1]);
          int i_node = node_coords[0]*(n_div + 1) + node_coords[1];
          data(i_node, all) << pnt.X(), pnt.Y(), pnt.Z();
          real(i_node) = 1;
        }
      }
      data *= 1e-3; // convert to meters
      Tecplot_file::Structured_block zone(file, n_div + 1, format_str(100, "surface%u", i_surf), 2);
      zone.write(data.data(), real.data());
    }
  } else {
    Tecplot_file file(file_name, 2, {}, 0.);
    for (unsigned i_curve = 0; i_curve < curves.size(); ++i_curve) {
      auto& curve = curves[i_curve];
      double param_bounds [2] {curve->FirstParameter(), curve->LastParameter()};
      Mat<dyn, dyn> data(n_div + 1, 2);
      for (int i_node = 0; i_node <= n_div; ++i_node) {
        double interp = double(i_node)/n_div;
        double param = (1 - interp)*param_bounds[0] + interp*param_bounds[1];
        auto pnt = curve->Value(param);
        data(i_node, all) << pnt.X(), pnt.Y();
      }
      data *= 1e-3; // convert to meters
      Tecplot_file::Structured_block zone(file, n_div + 1, format_str(100, "curve%u", i_curve), 1);
      zone.write(data.data(), nullptr);
    }
  }
}

Mat<> Occt_geom::nearest_point(Mat<> point)
{
  HEXED_ASSERT(point.size() == nd, format_str(100, "`point` must be %iD", nd));
  Mat<> scaled = point*1000; // convert to mm
  math::Nearest_point<> nearest(scaled);
  if (nd == 2) {
    gp_Pnt2d occt_point(scaled(0), scaled(1));
    // iterate through curves and find which ones has the nearest point
    for (auto& curve : curves) {
      Geom2dAPI_ProjectPointOnCurve proj(occt_point, curve);
      if(proj.NbPoints()) {
        gp_Pnt2d occt_candidate = proj.NearestPoint();
        nearest.merge(Mat<2>{occt_candidate.X(), occt_candidate.Y()});
      }
    }
  } else {
    gp_Pnt occt_point(scaled(0), scaled(1), scaled(2));
    // iterate through the surfaces and find which one has the nearest point
    for (auto& surface : surfaces) {
      GeomAPI_ProjectPointOnSurf proj(occt_point, surface);
      if (proj.IsDone()) {
        gp_Pnt occt_candidate = proj.NearestPoint();
        nearest.merge(Mat<3>{occt_candidate.X(), occt_candidate.Y(), occt_candidate.Z()});
      }
    }
  }
  return nearest.point()/1000;
}

std::vector<double> Occt_geom::intersections(Mat<> point0, Mat<> point1)
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
    // iterate through curves and compute the indersections with each
    for (auto& curve : curves) {
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
    // iterate through surfaces and compute the indersections with each
    for (auto& surface : surfaces) {
      // compute intersections
      GeomAPI_IntCS inter(curve, surface);
      // translate to our parametric format
      int n = inter.NbPoints();
      for (int i = 0; i < n; ++i) {
        double params [3];
        inter.Parameters(i + 1, params[0], params[1], params[2]);
        sects.push_back(params[2]/dist);
      }
    }
  }
  return sects;
}

void Occt_geom::write_image(const TopoDS_Shape& shape, std::string file_name, Mat<3> eye_pos, Mat<3> look_at_pos, int resolution)
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
TopoDS_Shape Occt_geom::execute_reader(std::string file_name)
{
  set_message();
  reader_t reader;
  auto result = reader.ReadFile(file_name.c_str());
  HEXED_ASSERT(result == IFSelect_RetDone, "failed to read geometry file");
  reader.TransferRoots();
  return reader.OneShape();
}

TopoDS_Shape Occt_geom::read(std::string file_name)
{
  unsigned extension_start = file_name.find_last_of(".");
  HEXED_ASSERT(extension_start != std::string::npos, "`file_name` has no extension");
  std::string case_sensitive = file_name.substr(extension_start + 1, std::string::npos);
  std::string ext = case_sensitive;
  for (char& c : ext) c = tolower(c);
  if      (ext == "igs" || ext == "iges") return execute_reader<IGESControl_Reader>(file_name);
  else if (ext == "stp" || ext == "step") return execute_reader<STEPControl_Reader>(file_name);
  throw std::runtime_error(format_str(1000, "`hexed::Occt_geom::read` failed to recognize file exteinsion `.%s`.", case_sensitive.c_str()));
}

std::vector<Mat<3, 3>> triangles(opencascade::handle<Poly_Triangulation> poly)
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

opencascade::handle<Poly_Triangulation> Occt_geom::read_stl(std::string file_name, double scale)
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

std::vector<Mat<2, 2>> segments(const TopoDS_Shape& shape, int n_segments)
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
