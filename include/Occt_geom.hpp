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
  // whether the message settings for OCCT have been set to prevent it polluting `std::cout`
  static bool message_set; // defaults to false
  // set OCCT messages to do nothing, if that hasn't already been done
  // this should be called at the start of any function that does any file IO
  static void set_message();
  TopoDS_Shape topo_shape;
  public:
  /*! \brief Construct directly from an OCCT shape object.
   * \details Coordinates are interpreted dimensionally.
   * OCCT works exclusively in mm and hexed works exclusively in m,
   * so the numerical value of the input shape's dimensions will be scaled by 1/1000.
   * That said, if you're using `read` to import a geometry from a file,
   * all the unit conversions happen automatically so you don't have to worry about it.
   */
  Occt_geom(TopoDS_Shape&&);
  Mat<> nearest_point(Mat<> point) override;
  //! \note May return duplicate points if intersection is on the boundary of multiple faces.
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override;
  inline TopoDS_Shape shape() {return topo_shape;} //!< fetches the underlying OCCT geometry representation
  /*! \brief Renders an image of the geometry and writes it to an image file.
   * \details For verifying/debuggin CAD translations.
   * \param file_name should include a file type extension which determines
   * the format of the image file.
   * I know that `.png` is supported, but I'm not sure what else,
   * and to be totally honest I don't really care.
   * If you're overcome by curiousity, feel free to go poke your nose around
   * `Image_AlienPixMap::Save` in OCCT.
   * \param resolution Width/height of image in pixels (it's always square).
   */
  void write_image(std::string file_name, std::vector<Mat<3>> rotations = {}, int resolution = 1000);

  /*! \brief Reads a CAD file and constructs an `Occt_geom` from it.
   * \details Supply an OCCT CAD reader type as the template argument.
   * Each reader can read exactly one file type,
   * so be sure the type of the file matches the reader you pick.
   * Currently supported readers are:
   * - `IGESControl_Reader` : IGES files
   *
   * Coordinates are interpreted dimensionally and automatically converted to m.
   * Not thread safe.
   */
  template <typename reader_t>
  static Occt_geom read(std::string file_name)
  {
    set_message();
    reader_t reader;
    auto result = reader.ReadFile(file_name.c_str());
    HEXED_ASSERT(result == IFSelect_RetDone, "failed to read geometry file");
    reader.TransferRoots();
    return Occt_geom(reader.OneShape());
  }
};

}
#endif
#endif
