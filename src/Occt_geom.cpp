#include <Occt_geom.hpp>
#if HEXED_USE_OCCT

// geometry
#include <TopoDS_Iterator.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS_Face.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GC_MakeLine.hxx>
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

namespace hexed
{

bool Occt_geom::message_set = false;

void Occt_geom::set_message()
{
  if (!message_set) {
    message_set = true;
    auto& printers = Message::DefaultMessenger()->ChangePrinters();
    printers.Clear();
    Handle(Message_PrinterToReport) printer;
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

Occt_geom::Occt_geom(const TopoDS_Shape& shape, int n_dim)
: nd{n_dim}
{
  // collect all surfaces
  if (nd == 3) {
    iterate(shape, TopAbs_FACE, [&](const TopoDS_Shape& s){
      TopoDS_Face face = TopoDS::Face(s);
      surfaces.push_back(BRep_Tool::Surface(face));
    });
  }
}

Mat<> Occt_geom::nearest_point(Mat<> point)
{
  HEXED_ASSERT(point.size() == nd, format_str(100, "`point` must be %iD", nd));
  Mat<> scaled = point*1000; // convert to mm
  // initialize nearest point and distance
  Mat<> nearest = scaled;
  double dist = std::numeric_limits<double>::max();
  if (nd == 3) {
    gp_Pnt occt_point(scaled(0), scaled(1), scaled(2));
    // iterate through the surfaces
    // and find which one has the nearest point
    for (auto& surface : surfaces) {
      // find the nearest point on this surface
      GeomAPI_ProjectPointOnSurf proj(occt_point, surface);
      gp_Pnt occt_candidate = proj.NearestPoint();
      Mat<3> candidate{occt_candidate.X(), occt_candidate.Y(), occt_candidate.Z()};
      // if this point is closer than the nearest point so far, make it the new nearest point
      double d = (candidate - scaled).norm();
      if (d < dist) {
        dist = d;
        nearest = candidate;
      }
    }
  }
  return nearest/1000;
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
  if (nd == 3) {
    // compute the line through the given points
    gp_Pnt pnt0(scaled0(0), scaled0(1), scaled0(2));
    gp_Pnt pnt1(scaled1(0), scaled1(1), scaled1(2));
    Handle(Geom_Line) line = GC_MakeLine(pnt0, pnt1);
    Handle(Geom_Curve) curve = line;
    // iterate through the surfaces and compute the indersections with each
    for (auto& surface : surfaces) {
      // compute intersections
      GeomAPI_IntCS inter(curve, surface);
      // translate to our parametric format
      int n = inter.NbPoints();
      for (int i = 0; i < n; ++i) {
        double params [3];
        inter.Parameters(i+1, params[0], params[1], params[2]);
        sects.push_back(params[2]/dist);
      }
    }
  }
  return sects;
}

void Occt_geom::write_image(const TopoDS_Shape& shape, std::string file_name, Mat<3> eye_pos, Mat<3> look_at_pos, int resolution)
{
  // general setup
  Handle(Aspect_DisplayConnection) displayConnection = new Aspect_DisplayConnection();
  Handle(OpenGl_GraphicDriver) graphicDriver = new OpenGl_GraphicDriver(displayConnection);
  Handle(V3d_Viewer) viewer = new V3d_Viewer(graphicDriver);
  viewer->SetDefaultLights();
  viewer->SetLightOn();
  Handle(AIS_InteractiveContext) context = new AIS_InteractiveContext(viewer);
  Handle(V3d_View) view = viewer->CreateView();
  Handle(Xw_Window) win = new Xw_Window(graphicDriver->GetDisplayConnection(), "", 0, 0, resolution, resolution);
  win->SetVirtual(true);
  view->SetWindow(win);
  view->SetBackgroundColor(Quantity_Color(Quantity_NOC_BLACK));
  view->MustBeResized();
  view->AutoZFit();
  // add shapes
  Handle(AIS_Shape) shaded = new AIS_Shape(shape);
  context->Display(shaded, false);
  context->SetDisplayMode(shaded, AIS_Shaded, false);
  Handle(AIS_Shape) wireframe = new AIS_Shape(shape);
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

}
#endif
