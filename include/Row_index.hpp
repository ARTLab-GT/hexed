#ifndef HEXED_ROW_INDEX_HPP_
#define HEXED_ROW_INDEX_HPP_

#include "math.hpp"

namespace hexed
{

/*! \brief Facilitates the process of iterating through all rows along a particular dimension.
 * \details This is essential in the local kernel where a `row_size`[x`row_size`[x`row_size`]]
 * array of flux values must be multipled by the derivative matrix row-wise
 * along every dimension.
 */
class Row_index
{
  int i_outer = 0;
  int i_inner = 0;
  int i_fq = 0;

  public:
  const int n_dim;
  const int row_size;
  const int i_dim;
  const int n_qpoint; //!< \brief how many total quadrature points an element has
  const int n_fqpoint; //!< \brief how many quadrature points the face of an element has
  const int stride; //!< \brief the stride required to step through the row

  /*! constructs a `Row_index` in `nd` dimensions, with row size `rs`,
   * and applicable to rows in the `id`th dimension.
   * Object is initialized to point to the first row.
   */
  constexpr Row_index(int nd, int rs, int id)
  : n_dim{nd}, row_size{rs}, i_dim{id},
    n_qpoint{math::pow(row_size, n_dim)},
    n_fqpoint{n_qpoint/row_size},
    stride{math::pow(row_size, n_dim - 1 - i_dim)}
  {}

  //! \brief advances object to the next row
  constexpr void operator++()
  {
    ++i_fq;
    ++i_inner;
    if (i_inner == stride) {
      i_inner = 0;
      ++i_outer;
    }
  }

  //! \brief returns `true` if object is pointing to a valid row
  //! and `false` if it is pointing past the last row
  constexpr operator bool() const {return i_fq < n_fqpoint;}

  //! \brief returns the index of the row (or equivalently the face quadrature point) the object is pointing to
  constexpr int i_face_qpoint() const {return i_fq;}

  //! \brief returns the index (with respect to the entire quadrature point array) of
  //! the `i_node`th quadrature point in the current row
  constexpr int i_qpoint(int i_node) const {return i_outer*stride*row_size + i_inner + i_node*stride;}

  //! \brief returns the index of a quadrature point in its own row, i.e. which `Basis::node` it is
  //! \details always gives the same value regardless of which `i_face_qpoint()` this object is pointing to
  constexpr int i_node(int i_qpoint) const {return (i_qpoint/stride)%row_size;}
};

}
#endif
