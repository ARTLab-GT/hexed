#ifndef HEXED_ROW_INDEX_HPP_
#define HEXED_ROW_INDEX_HPP_

#include "math.hpp"

namespace hexed
{

class Row_index
{
  int i_outer = 0;
  int i_inner = 0;
  int i_fq = 0;

  public:
  const int n_dim;
  const int row_size;
  const int i_dim;
  const int n_qpoint;
  const int n_fqpoint;
  const int stride;

  constexpr Row_index(int nd, int rs, int id)
  : n_dim{nd}, row_size{rs}, i_dim{id},
    n_qpoint{custom_math::pow(row_size, n_dim)},
    n_fqpoint{n_qpoint/row_size},
    stride{custom_math::pow(row_size, n_dim - 1 - i_dim)}
  {}

  constexpr void operator++()
  {
    ++i_fq;
    ++i_inner;
    if (i_inner == stride) {
      i_inner = 0;
      ++i_outer;
    }
  }

  constexpr operator bool() const {return i_fq < n_fqpoint;}
  constexpr int i_face_qpoint() const {return i_fq;}
  constexpr int i_qpoint(int i_node) const {return i_outer*stride*row_size + i_inner + i_node*stride;}
};

}
#endif
