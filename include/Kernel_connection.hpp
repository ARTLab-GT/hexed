#ifndef HEXED_KERNEL_CONNECTION_HPP_
#define HEXED_KERNEL_CONNECTION_HPP_

#include <utils.hpp>

namespace hexed {

class Connection_direction {
  public:
  std::array<int, 2> i_dim;
  std::array<bool, 2> face_sign;
  int rotate = 0;
  inline int i_face(int i_side) const {return 2*i_dim[i_side] + face_sign[i_side];}
  /*!
   * Answers the question: Is it necessary to flip the normal of element `i_side` so that it
   * points from element 0 into element 1?
   */
  inline bool flip_normal(int i_side) const {return face_sign[i_side] == i_side;}
  /*!
   * Answers the question: Is it neccesary to flip axis `face_index(0).i_dim` of element 1
   * to match the coordinate systems?
   */
  inline bool flip_tangential() const {
    //! if you're swapping two axes, you have to flip one of them to make a valid rotation. If you're not
    //! flipping a normal (or flipping both of them) then you have to flip a tangential
    return (i_dim[0] != i_dim[1]) && (flip_normal(0) == flip_normal(1));
  }
  /*!
   * Answers the question: Is it necessary to transpose the rows/columns of the face
   * quadrature points of element 1 to match element 0? Only applicable to 3D, where some
   * face combinations can create a row vs column major mismatch. If 2D, always returns `false`.
   */
  inline bool transpose() const {
    return ((i_dim[0] == 0) && (i_dim[1] == 2)) || ((i_dim[0] == 2) && (i_dim[1] == 0));
  }
};

//! \relates Connection_direction
bool operator==(Connection_direction dir0, Connection_direction dir1);
//! \relates Connection_direction
std::string to_string(Connection_direction);

struct Hard_kernel_connection {
  Connection_direction direction;
  double* state [2][2];
  double* normal;
  int mask [2];
  double nominal_area;
};

struct Kernel_face_refinement {
  public:
  double* coarse [2];
  double* fine [2][2];
  int mask;
  int split_dim;
};

}
#endif
