#ifndef HEXED_KERNEL_CONNECTION_HPP_
#define HEXED_KERNEL_CONNECTION_HPP_

namespace hexed
{

class Connection_direction
{
  public:
  std::array<int, 2> i_dim;
  std::array<bool, 2> face_sign;
  int i_face(int i_side) {return 2*i_dim[i_side] + face_sign[i_side];}
  /*!
   * Answers the question: Is it necessary to flip the normal of element `i_side` so that it
   * points from element 0 into element 1?
   */
  bool flip_normal(int i_side) {return face_sign[i_side] == i_side;}
  /*!
   * Answers the question: Is it neccesary to flip axis `face_index(0).i_dim` of element 1
   * to match the coordinate systems?
   */
  bool flip_tangential()
  {
    //! if you're swapping two axes, you have to flip one of them to make a valid rotation. If you're not
    //! flipping a normal (or flipping both of them) then you have to flip a tangential
    return (i_dim[0] != i_dim[1]) && (flip_normal(0) == flip_normal(1));
  }
  /*!
   * Answers the question: Is it necessary to transpose the rows/columns of the face
   * quadrature points of element 1 to match element 0? Only applicable to 3D, where some
   * face combinations can create a row vs column major mismatch. If 2D, always returns `false`.
   */
  bool transpose()
  {
    return ((i_dim[0] == 0) && (i_dim[1] == 2)) || ((i_dim[0] == 2) && (i_dim[1] == 0));
  }
};

//! \brief Represents a connection between elements as the kernel sees it.
//! \details Similar idea to `Kernel_element`.
class Kernel_connection
{
  public:
  virtual Connection_direction get_direction() = 0;
  virtual double* state() = 0; //!< \brief state data for both sides \details layout: [is_ldg_state][i_side][i_var][i_face_qpoint]
  //! \brief state data for one side and for either the extrapolated state or the LDG storage
  //! \details layout: [i_var][i_face_qpoint]
  virtual double* new_state(int i_side, bool is_ldg) = 0;
  virtual double* normal() = 0; //!< \brief face normal vector \details `nullptr` for Cartesian
};

}
#endif
