#ifndef HEXED_VERTEX_HPP_
#define HEXED_VERTEX_HPP_

#include <array>
#include <vector>
#include <memory>
#include <set>
#include <unordered_set>
#include <Eigen/Dense>
#include "Lock.hpp"

namespace hexed
{

/*! \brief Represents a vertex in a deformed grid.
 * \details Used by `Deformed_element`s to manage data
 * that is shared by nodal neighbors, including position and artificial viscosity coefficient.
 * \note Construction and ownership should be managed by the `Transferable_ptr` member class.
 */
class Vertex
{
  public:
  class Transferable_ptr;
  class Non_transferable_ptr;
  typedef double (Eigen::VectorXd::*reduction)() const;
  static constexpr reduction vector_max = &Eigen::VectorXd::maxCoeff;
  static constexpr reduction vector_min = &Eigen::VectorXd::minCoeff;
  std::array<double, 3> pos {0, 0, 0};
  bool mobile = false;
  std::vector<int> record; //!< for algorithms to keep notes as they please
  Lock lock; //!< for any algorithms that could involve data races on vertices

  ~Vertex();
  //! if we have a reason to copy/move vertices, we can implement these
  Vertex(const Vertex&) = delete;
  Vertex(Vertex&&) = delete;
  Vertex& operator=(const Vertex&) = delete;
  Vertex& operator=(Vertex&&) = delete;
  int mass();
  /*! \brief Specify that another vertex represents the same grid point as `*this`.
   * \details `*this` will acquire `other`'s resources:
   * - All `Transferable_ptr`s pointing to `other` will be changed to point to `this`.
   * - All `Non_transferable_ptr`s pointing to `other` will be `nullify`d.
   * - The mass of `other` will be added to that of `*this`.
   * - The `pos`s will be averaged, weighted by `mass`.
   * - The vertex will be mobile if both vertices were mobile.
   * The fact that "eat" seemed like the obvious word to describe this might be a
   * sign that I've been reading too much SnK...
   */
  void eat(Vertex& other);
  void calc_relax(double factor = .5); //!< compute a (but do not apply) new position resulting in a smoother grid.
  void apply_relax(); //!< update `pos` to the position computed by `calc_relax`.
  //! determine the shared shareable_value of all `Transferable_ptr`s to this by applying `reduction`. Thread safe.
  double shared_value(reduction = vector_max); //!< cppcheck-suppress internalAstError
  static void connect(Vertex&, Vertex&); //!< specify that two vertices are connected by an edge
  static bool are_neighbors(Vertex&, Vertex&);

  private:
  int m;
  std::array<double, 3> relax {0, 0, 0};
  std::unordered_set<Transferable_ptr*> trbl_ptrs;
  std::unordered_set<Non_transferable_ptr*> nont_ptrs;
  std::set<Vertex*> neighbors; // use set because have to iterate through
  Vertex(std::array<double, 3> pos);
};

/*! \brief A pointer class which can be transferred by `Vertex::eat`.
 * \details This is used by
 * `Deformed_element` to track its vertices, which are initially separate and then
 * eat those of neighboring elements as the grid connectivity is established.
 */
class Vertex::Transferable_ptr
{
  std::shared_ptr<Vertex> ptr;

  public:
  //! use this member to store a value that should be synchronized among all elements sharing this vertex
  double shareable_value;

  //! Create a vertex and construct a `Transferable_ptr` to it.
  Transferable_ptr(std::array<double, 3> pos);
  //! copy semantics creates another pointer to the same vertex
  Transferable_ptr(const Transferable_ptr&);
  ~Transferable_ptr();
  Transferable_ptr& operator=(const Transferable_ptr&);

  Vertex* operator->();
  const Vertex* operator->() const;
  Vertex& operator*();
  const Vertex& operator*() const;
};

/*! \brief A pointer class which is not transferred by `Vertex::eat`.
 * \details This is used by `Deformed_grid` to maintain a list of vertices without duplicates.
 */
class Vertex::Non_transferable_ptr
{
  Vertex* ptr;

  public:
  Non_transferable_ptr(Vertex&);
  Non_transferable_ptr(const Non_transferable_ptr&);
  ~Non_transferable_ptr();
  Non_transferable_ptr& operator=(const Non_transferable_ptr&);

  //! return `true` if vertex `this` points to exists and `false` if it does not
  operator bool() const;
  //! Warning: if `operator bool` is `false` the following will return `nullptr`!
  Vertex* operator->();
  Vertex& operator*();
  //! assert that the vertex `this` points to no longer exists. `operator bool` will now return `false`
  void nullify();
  static inline bool is_null(const Non_transferable_ptr& ntptr) {return !ntptr;}
};

}

#endif
