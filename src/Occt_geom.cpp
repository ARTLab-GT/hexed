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
  TopoDS_Iterator it(topo_shape);
  Mat<3> nearest = scaled;
  double dist = std::numeric_limits<double>::max();
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
  return {};
}

}
#endif
