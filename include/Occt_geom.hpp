#ifndef HEXED_OCCT_GEOM_HPP_
#define HEXED_OCCT_GEOM_HPP_

#include "config.hpp"
#if HEXED_USE_OCCT

#include <TopoDS.hxx>
#include <IGESControl_Reader.hxx>
#include "Surface_geom.hpp"

namespace hexed
{

/*! \brief Interface for geometry defined by
 * [Open CASCADE Technology](https://dev.opencascade.org/doc/overview/html/index.html) (OCCT).
 * \details Allows you to mesh a geometry defined by OCCT
 * by implementing the `Surface_geom` functionality on top of the OCCT API.
 * The main application for this is reading geometry from CAD files,
 * but this could also be used to pipe in CAD directly from another OCCT program.
 * \warning For now, only supports 3D, so all input points must be 3D.
 * \todo implement 2D.
 */
class Occt_geom : public Surface_geom
{
  public:
  //! \see `read_cad()`
  Occt_geom(TopoDS_Shape&&);
  Mat<> nearest_point(Mat<> point) override;
  //! \note May return duplicate points if intersection is on the boundary of multiple faces.
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override;
};

/*! \brief Reads a CAD file and returns a geometry object.
 * \details Returned object can be used to construct an `Occt_geom`.
 * Supply an OCCT CAD reader type as the template argument.
 * Each reader can read exactly one file type,
 * so be sure the type of the file matches the reader you pick.
 * Currently supported readers are:
 * - `IGESControl_Reader` : IGES files
 */
template <typename reader_t>
TopoDS_Shape read_cad(std::string file_name)
{
  reader_t reader;
  auto result = reader.ReadFile(file_name.c_str());
  HEXED_ASSERT(result == IFSelect_RetDone, "failed to read geometry file");
  reader.TransferRoots();
  return reader.OneShape();
}

}
#endif
#endif
