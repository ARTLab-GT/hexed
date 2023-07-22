#ifndef HEXED_MESH_INTERFACE_HPP_
#define HEXED_MESH_INTERFACE_HPP_

#include <string>
#include "criteria.hpp"
#include "Boundary_condition.hpp"
#include "Layer_sequence.hpp"
#include "Surface_geom.hpp"

namespace hexed
{

class Mesh_interface
{
  public:
  //! \returns Nominal size (\f$\Delta h\f$) of elements with refinement level 0.
  virtual double root_size() = 0;

  //! \name Automated tree meshing
  //! \{
  /*! \brief Initializes meshing with bin/quad/[octree](https://en.wikipedia.org/wiki/Octree) topology.
   * \details Creates a tree initialized with one element with ref level 0.
   * Technically, free-form elements can also be created in the same mesh,
   * but they will not be connected to the tree.
   * Only one tree can be created.
   * Connections and boundary conditions are set automatically for tree elements,
   * so tree meshes should always automatically be valid unless you explicitly invalidate it with `Mesh::disconnect_boundary`.
   * \param extremal_bcs Boundary conditions to apply to elements
   *   with exposed faces at the extremal boundaries of the tree.
   *   Takes ownership of objects that the vector entries point to.
   *   It must contain exactly `2*n_dim` entries.
   *   E.g. `extremal_bcs[0]` is the boundary condition to apply to the minimum \f$x_0\f$ face,
   *   `extremal_bcs[1]` is the BC for the maximum \f$x_0\f$ face,
   *   `extremal_bcs[2]` is for minimum \f$x_1\f$.
   * \param origin The minimal corner of the tree root will be located at `origin`.
   */
  virtual void add_tree(std::vector<Flow_bc*> extremal_bcs, Mat<> origin = Mat<>::Zero(3)) = 0;
  /*! \brief Defines the surface geometry to be meshed as a boundary.
   * \details Acquires ownership of the objects pointed to `geometry` and `surface_bc`.
   * The geometric surface represented by `geometry` is now a boundary of the domain
   * and its boundary condition is `surface_bc`.
   * The flood fill algorithm is then executed starting at `flood_fill_start`
   * to determine which elements are in the domain
   * and extrusion is executed to create a suitably body-fitted mesh.
   * `flood_fill_start` must have at least `n_dim` entries and extra entries are ignored.
   * If the mesh is excessively coarse, there may be no elements in the domain
   * as they are all too close to the surfaces.
   * So, `Mesh::update` should be calling a few times before `Mesh::set_surfaces`.
   * Any surfaces defined by previous invokations of `set_surface` are forgotten.
   */
  virtual void set_surface(Surface_geom* geometry, Flow_bc* surface_bc, Eigen::VectorXd flood_fill_start = Eigen::VectorXd::Zero(3)) = 0;
  /*! \brief Updates tree mesh based on user-supplied (un)refinement criteria.
   * \details Evaluates `refine_criterion` and `unrefine_criterion` on every element in the tree (if a tree exists).
   * Whenever `refine_criterion` is `true` and `unrefine_criterion` is `false`, that element is refined.
   * Whenever `unrefine_criterion` is `true` and `refine_criterion` is `false` for a complete group of sibling elements,
   * that group is unrefined.
   * In order to satisfy some criteria regarding the refinement level of neighbors,
   * some additional elements may be refined and some elements may not be unrefined.
   * Both criteria must be thread-safe
   * and must not depend on the order in which elements are processed.
   * The flood fill and extrusion are also updated.
   * \returns `true` if the mesh was changed, else `false`
   */
  virtual bool update(std::function<bool(Element&)> refine_criterion = criteria::always, std::function<bool(Element&)> unrefine_criterion = criteria::never) = 0;
  //! \brief Relax the vertices to improve mesh quality.
  //! \param factor A larger number yields more change in the mesh. 0 => no update, 1 => "full" update, > 1 allowed but suspect
  virtual void relax(double factor = 0.9) = 0;
  virtual int surface_bc_sn() = 0; //!< what is the serial number of the geometry surface BC?
  //! \}

  //! \name observers
  //!\{
  virtual int n_elements() = 0; //!< \brief number of elements currently in the mesh
  //! An object to provide information about whether the mesh connectivity is valid and if not, why.
  class Connection_validity
  {
    public:
    const int n_redundant; //!< number of redundant connections, counting each participating face as one
    const int n_missing; //!< number of missing connections
    //! returns true if connectivity is valid
    inline operator bool() {return (n_redundant == 0) && (n_missing == 0);}
    //! if connectivity is invalid, throw an exception with a helpful message
    inline void assert_valid()
    {
      if (!*this) {
        auto message = "Invalid mesh with " + std::to_string(n_redundant) + " redundant connections and "
                       + std::to_string(n_missing) + " missing connections.";
        throw std::runtime_error(message);
      }
    }
  };
  /*! \brief Returns a `Connection_validity` object describing whether the mesh connectivity is valid.
   * \details Suggested uses:
   * - `if (mesh.valid()) {\\...do something that requires a valid mesh}`
   * - `mesh.valid().assert_valid();`
   */
  virtual Connection_validity valid() = 0;
  //! minimum information required to identify an element
  struct elem_handle {int ref_level; bool is_deformed; int serial_n;};
  //! get handles for all elements currently in the mesh, in no particular order (mostly for testing/debugging)
  virtual std::vector<elem_handle> elem_handles() = 0;
  //! Temporarily resets the vertices of a mesh to their nominal positions for debugging
  class Reset_vertices
  {
    Mesh_interface& m;
    public:
    //! resets to nominal position
    inline Reset_vertices(Mesh_interface& mesh) : m{mesh} {m.reset_verts();}
    //! restores vertices to where they were before this object was constructed
    inline ~Reset_vertices() {m.restore_verts();}
  };
  //!\}

  protected:
  virtual void reset_verts() = 0;
  virtual void restore_verts() = 0;
};

}
#endif
