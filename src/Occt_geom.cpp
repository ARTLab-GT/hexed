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

Occt_geom::Occt_geom(TopoDS_Shape&& s)
: topo_shape{s}
{}

Mat<> Occt_geom::nearest_point(Mat<> point)
{
  HEXED_ASSERT(point.size() == 3, "`point` must be 3D");
  Mat<3> scaled = point*1000;
  gp_Pnt occt_point(scaled(0), scaled(1), scaled(2));
  Mat<3> nearest = scaled;
  double dist = std::numeric_limits<double>::max();
  TopoDS_Iterator it(topo_shape);
  while (it.More()) {
    if (it.Value().ShapeType() == TopAbs_FACE) {
      TopoDS_Face face = TopoDS::Face(it.Value());
      Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
      GeomAPI_ProjectPointOnSurf proj(occt_point, surface);
      gp_Pnt occt_candidate = proj.NearestPoint();
      Mat<3> candidate{occt_candidate.X(), occt_candidate.Y(), occt_candidate.Z()};
      double d = (candidate - scaled).norm();
      if (d < dist) {
        dist = d;
        nearest = candidate;
      }
    }
    it.Next();
  }
  return nearest/1000;
}

std::vector<double> Occt_geom::intersections(Mat<> point0, Mat<> point1)
{
  HEXED_ASSERT(point0.size() == 3, "`point0` must be 3D");
  HEXED_ASSERT(point1.size() == 3, "`point1` must be 3D");
  Mat<3> scaled0 = 1000*point0;
  Mat<3> scaled1 = 1000*point1;
  double dist = (scaled1 - scaled0).norm();
  gp_Pnt pnt0(scaled0(0), scaled0(1), scaled0(2));
  gp_Pnt pnt1(scaled1(0), scaled1(1), scaled1(2));
  Handle(Geom_Line) line = GC_MakeLine(pnt0, pnt1);
  Handle(Geom_Curve) curve = line;
  std::vector<double> sects;
  TopoDS_Iterator it(topo_shape);
  while (it.More()) {
    if (it.Value().ShapeType() == TopAbs_FACE) {
      TopoDS_Face face = TopoDS::Face(it.Value());
      Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
      GeomAPI_IntCS inter(curve, surface);
      int n = inter.NbPoints();
      for (int i = 0; i < n; ++i) {
        double params [3];
        inter.Parameters(i+1, params[0], params[1], params[2]);
        sects.push_back(params[2]/dist);
      }
    }
    it.Next();
  }
  return sects;
}

void Occt_geom::write_image(std::string file_name, std::vector<Mat<3>> rotations, int resolution)
{
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
  Handle(AIS_Shape) presentation = new AIS_Shape(topo_shape);
  context->Display(presentation, Standard_False);
  context->SetDisplayMode(presentation, AIS_Shaded, Standard_False);
  view->SetProj(V3d_Zneg);
  for (Mat<3> rot : rotations) {
    view->Rotate(rot(0), rot(1), rot(2), 0, 0, 0, 0);
  }
  view->FitAll(.2);
  view->Redraw();
  view->Dump(file_name.c_str());
}

}
#endif
