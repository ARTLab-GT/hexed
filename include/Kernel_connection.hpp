#ifndef HEXED_KERNEL_CONNECTION_HPP_
#define HEXED_KERNEL_CONNECTION_HPP_

#include <utils.hpp>

namespace hexed {

//! \brief Describes the direction in which to faces are connected.
//! \details In other words, it describes the relative orientation of the connected faces in reference coordinates
//! so that a consistent mapping between their quadrature points and normals can be established.
class Connection_direction {
  public:
  std::array<int, 2> i_dim; //!< \brief The normal axis of each face in reference coordinates
  //! \brief The sign of the nonzero component of the normal vector of each face in reference coordinates
  std::array<bool, 2> face_sign;
  //! \brief Indicates whether the coordinates need to be rotated to make the quadrature points match.
  //! \details Specifically, int indicates that face 1 is rotated by `rotate` 90\f$\degree\f$ increments
  //! about the shared normal relative to face 0.
  //! So, when permuting quadrature points or vertices, after all other transformations have been performed,
  //! `rotate` 90\f$\degree\f$ rotations will be performed on face 1 __in a negative sense (i.e. clockwise)__
  //! to undo the specified rotation.
  int rotate = 0;
  inline int i_face(int i_side) const {return 2*i_dim[i_side] + face_sign[i_side];}
  //! \brief Returns `true` if the normal of face `i_side` must be inverted to obey the sign convention.
  //! \details The sign convention is for the normal to point from the element on side 0 into that on side 1.
  inline bool flip_normal(int i_side) const {return face_sign[i_side] == i_side;}
  //! \brief Returns `true` if the indices of axis 0 of face 1 must be reversed to match the indexing of face 0.
  inline bool flip_tangential() const {
    // if you're swapping two axes, you have to flip one of them to make a valid rotation. If you're not
    // flipping a normal (or flipping both of them) then you have to flip a tangential
    return (i_dim[0] != i_dim[1]) && (flip_normal(0) == flip_normal(1));
  }
  //! \brief Returns `true` if the rows and columns of face 1 must be transposed to match the indexing of face 0.
  //! \details In 2D, this is always `false`.
  inline bool transpose() const {
    return ((i_dim[0] == 0) && (i_dim[1] == 2)) || ((i_dim[0] == 2) && (i_dim[1] == 0));
  }
};

//! \relates Connection_direction
bool operator==(Connection_direction dir0, Connection_direction dir1);
//! \relates Connection_direction
std::string to_string(Connection_direction);
//! \relates Connection_direction
std::ostream& operator<<(std::ostream&, Connection_direction);

//! \brief The minimal amount of data that the kernels need about a `Neighbor_connection`
struct Hard_kernel_connection {
  Connection_direction direction;
  double* state [2][2];
  double* normal;
  int mask [2];
  double nominal_area;
};

//! \brief The minimal amount of data that the kernels need about `Face_refinement`
struct Kernel_face_refinement {
  public:
  //! \brief Pointers to the convective and ldg storage arrays for the coarse face. Layout: [is_ldg]
  double* coarse [2];
  //! \brief Pointers to the convective and ldg storage arrays for the fine faces. Layout: [i_side][is_ldg]
  double* fine [2][2];
  int mask; //! \brief The mask level of the faces (which is the same)
  int split_dim; //! \brief Whether the face is split along dimension 0 or (in 3D) dimension 1
};

}
#endif
