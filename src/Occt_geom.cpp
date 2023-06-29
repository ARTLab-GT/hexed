#include <Occt_geom.hpp>
#if HEXED_USE_OCCT

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

Occt_geom::Occt_geom(TopoDS_Shape&&)
{}

Mat<> Occt_geom::nearest_point(Mat<> point)
{
  return Mat<>::Zero(3);
}

std::vector<double> Occt_geom::intersections(Mat<> point0, Mat<> point1)
{
  return {};
}

}
#endif
