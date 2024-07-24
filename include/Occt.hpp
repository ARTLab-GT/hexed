#ifndef HEXED_OCCT_HPP_
#define HEXED_OCCT_HPP_

#include "config.hpp"
#if HEXED_USE_OCCT

#include <TopoDS.hxx>
#include <Geom_Surface.hxx>
#include <Geom2d_Curve.hxx>
#include <Poly_Triangulation.hxx>
#include "constants.hpp"
#include "Simplex_geom.hpp"

namespace hexed {

/*! \brief Wrapper interface for CAD geometry defined by
 * [Open CASCADE Technology](https://dev.opencascade.org/doc/overview/html/index.html) (OCCT).
 * \details The idea is to make it straightforward to take OCCT objects
 * and turn them into objects that Hexed can work with
 * without needing to understand the complexities of what they really mean in OCCT thinking.
 * The main application for this is reading geometry from CAD files,
 * but this could also be used to pipe in CAD directly from another program that uses OCCT.
 * For usage purposes, this class is just a namespace.
 * Don't actually instantiate it.
 * (In fact, you shouldn't be able to, since the constructor is deleted).
 * It is a class and not a namespace because it has private static data members
 * to facilitate startup tasks for interacting with OCCT.
 *
 * \section occt_units Units
 * All length/position values are interpreted dimensionally.
 * OCCT works exclusively in mm and hexed \ref units "works exclusively in m",
 * so when converting between hexed and OCCT objects all numerical values will be scaled by \f$ 10^{\pm3} \f$.
 * That said, unless you're manually interacting with the dimensions/coordinates of the OCCT objects,
 * all the unit conversions are automatic so you don't have to worry about it.
 *
 * \warning This class has memory leaks! I am tempted to blame that on the OCCT developers,
 * who seem to be a little old-fasioned about their memory management,
 * but maybe I'm just using it wrong.
 * Why is CAD software always a mess?
 * Maybe engineering is just always a mess...
 */
class Occt {
  // whether the message settings for OCCT have been set to prevent it polluting `std::cout`
  static bool message_set; // defaults to false
  // set OCCT messages to do nothing, if that hasn't already been done
  // this should be called at the start of any function that does any file IO
  static void set_message();
  // reads a file of a specific type
  template<typename reader_t> static TopoDS_Shape execute_reader(std::string file_name);

  struct Triangulation {
    std::vector<Mat<3, 3>> tris;
    std::vector<Mat<2, 3>> params;
    std::vector<int> faces;
  };
  static Triangulation _triangulate(opencascade::handle<Poly_Triangulation>, int face);
  static Triangulation _triangulate(TopoDS_Shape shape, double angle = 10*constants::degree, double deflection = huge);

  public:
  Occt() = delete; //!< \brief Don't instantiate this class. Use its static members.

  /*! \brief Reads a CAD file.
   * \details
   * This can then be discretized with `triangles(TopoDS_Shape, double, double)`.
   * Coordinates are interpreted dimensionally and automatically converted.
   * Not thread safe.
   * File format is inferred from the file extension.
   * File extension is case insensitive.
   * Supported formats and extensions are:
   * - IGES: `.igs`, `.iges`
   * - STEP: `.stp`, `.step`
   * \warning I've gotten a memory leak of over **40 MB** when calling this on a nonexistent file path
   * (accoring to LeakSanitizer anyway), so check your paths in advance!
   */
  static TopoDS_Shape read(std::string file_name);

  /*! \brief Reads an STL file and returns a triangulation object.
   * \details File can be in ASCII or binary format.
   * STL files contain no unit information, so it will be assumed to be in meters by default.
   * It will then be converted to mm for consistency with other OCCT objects.
   * If the STL was written with a different unit in mind,
   * you can supply that unit definition from `constants.hpp` as the `scale` argument
   * and it will be effectively converted from that unit.
   * In general, the geometry coordinates will be multiplied by `scale`.
   * Result can be piped to `triangles()` to ultimately construct a `Simplex_geom<3>`.
   */
  static opencascade::handle<Poly_Triangulation> read_stl(std::string file_name, double scale = constants::meter);

  /*! \brief Creates an array of simplices that can be used to construct a `Simplex_geom<3>`.
   * \details A `Poly_Triangulation` can be obtained from, e.g. `read_stl()`.
   */
  static std::vector<Mat<3, 3>> triangles(opencascade::handle<Poly_Triangulation>);
  /*! \brief Obtains a triangulation of a CAD geometry.
   * \details This is now the preferred way to interact with CAD geometry---`Occt_geom` instances are unreliable.
   * Size of the mesh is determined by `angle` and `deflection`,
   * where in both cases a smaller value results in a finer mesh.
   * Usually, it is preferable to use only `angle`,
   * since `deflection` doesn't do as well at refining high-curvature regions.
   * The result can be piped to `triangles()` to fetch the elements of the triangulation.
   * \param shape CAD object to triangulate. Can be obtained from, e.g., `read()`.
   * \param angle Max angle allowed between neighboring triangles (note that as always, angles are radian).
   * \param deflection Max distance allowed between points on the triangulation and the true surface.
   *        Note that this is a dimensional parameter
   *        and appropriate values will depend on the length scale of your simulation,
   *        which is why the only reasonable default is an irrelevantly large parameter.
   */
  static std::vector<Mat<3, 3>> triangles(TopoDS_Shape shape, double angle = 10*constants::degree,
                                          double deflection = huge);

  static Simplex_geom<3> triangulate(TopoDS_Shape shape, double angle = 10*constants::degree,
                                     double deflection = huge, int n_div_edge = 1000);

  /*! \brief Discretizes the curves in a `TopoDS_Shape` into segments of a polygonal line.
   * \details A `TopoDS_Shape` can be obtained from `read()`,
   * and the results can be used to construct a `Simplex_geom<2>`.
   */
  static std::vector<Mat<2, 2>> segments(const TopoDS_Shape&, int n_segments);

  /*! \brief A `Surface_geom` that interacts with a CAD object directly.
   * \details Represents a CAD object defined with the OCCT interface as a `Surface_geom`.
   * Implements `nearest_point` and `intersections`
   * by directly working with the OCCT projection and intersection functions,
   * which can in principle be faster than working with triangulations when high accuracy is desired.
   * However, it has two important drawbacks:
   * - Projections and intersections are not very robust,
   *   which defeats the whole purpose of Hexed which is robust and automated meshing.
   * - It cannot understand face trimming, which is required to correctly represent some 3D shapes defined with
   *   [BRep Topology](https://dev.opencascade.org/doc/occt-6.7.0/overview/html/user_guides__modeling_data.html#occt_modat_5_2_1).
   *   It only supports curves and faces that are purely parametric (they are not trimmed by any curves).
   *
   * As a result, the preferred way to interact with CAD geometry is to discretize it with `triangles()`
   * and then convert it to a `Simplex geom`.
   */
  class Geom : public Surface_geom {
    const int nd;
    std::vector<opencascade::handle<Geom_Surface>> _surfaces;
    std::vector<opencascade::handle<Geom2d_Curve>> _curves;
    std::unique_ptr<Simplex_geom<2>> _simplex2;
    std::unique_ptr<Simplex_geom<3>> _simplex3;
    Simplex_geom_nd* _simplex;
    public:
    /*! \brief Construct directly from an OCCT shape object.
     * \details The shape is interpreted to have dimensionality specified by `n_dim`,
     * which may be either 2 or 3.
     * All input points must have `n_dim` elements, as will all output points.
     * If 3D, only faces are considered.
     * If 2D, only edges are considered, and all are projected onto the \f$ x_2 = 0 \f$ plane
     * (i.e. xy-plane).
     * Coordinates are interpreted dimensionally and automatically converted to m.
     */
    Geom(const TopoDS_Shape&, int n_dim, double angle = 10*constants::degree,
         double deflection = huge, int n_segments = 1000);
    Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
    //! \note May return duplicate points if intersection is on the boundary of multiple faces.
    std::vector<double> intersections(Mat<> point0, Mat<> point1) override;
    //! \brief Visualizes the triangulated geometry
    void visualize(std::string format, std::string file_name);
  };
};

}
#endif
#endif
