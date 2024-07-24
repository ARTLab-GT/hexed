#ifndef HEXED_VERTEX_INDS_HPP_
#define HEXED_VERTEX_INDS_HPP_

#include <vector>
#include <array>
#include "Kernel_connection.hpp"

namespace hexed
{

//! \brief The indices required to permute the vertices of face 1 of a connection to match face 0.
std::vector<int> face_vertex_inds(int n_dim, const Connection_direction&);

/*! \brief The indices of the element vertices which participate in a deformed connection.
 * \details Ordering is such that vertices which align in physical space correspond in the lists
 * and the vertices of face 0 are ordered in the same way as they would be if the face were considered in isolation.
 */
std::array<std::vector<int>, 2> vertex_inds(int n_dim, const Connection_direction&);

}
#endif
