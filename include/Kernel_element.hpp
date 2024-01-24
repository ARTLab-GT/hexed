#ifndef HEXED_KERNEL_ELEMENT_HPP_
#define HEXED_KERNEL_ELEMENT_HPP_

namespace hexed
{

/*! \brief Represents an element as the kernels see it.
 * \details Includes only the bare minimum of information needed to run the \ref Spatial "spatial discretization kernels".
 * This minimizes the number of times where an irrelevant change necessitates recompiling the kernels,
 * which is expensive and thus very annoying.
 * It also helps to better define what the kernels do and what exactly their inputs and outputs are.
 * The kernel is responsible for knowing what the `Storage_params` and storage order of the element are.
 */
class Kernel_element
{
  public:
  //! \brief pointer to the data where the state variables are stored
  //! \details includes any non-conservation variables such as artificial viscosity coefficient
  virtual double* state() = 0;
  virtual double* residual_cache() = 0; //! \brief where the convective residual is stored for 2-stage time integration
  virtual double* time_step_scale() = 0; //!< \brief storage for local time step \details layout: [i_qpoint]
  virtual double& vertex_time_step_scale(int i_vertex) = 0; //!< \todo the kernel should not need this
  virtual double nominal_size() = 0; //!< \brief nominal edge length of the element before any vertex motion
  virtual double* face(int i_face) = 0; //! \brief where the extrapolated face data is stored
  virtual bool deformed() const = 0; //! \brief whether this element is deformed
  /*!
   * \details the `j_dim`th component (in physical space) of the normal vector of the level surface
   * of the `i_dim`th reference coordinate which passes through the `i_qpoint`th quadrature point.
   * the magnitude of the normal vector is weighted by the surface area in physical space.
   * equivalent definition (note `i_dim`, `j_dim` transposed):
   * \f[ J^{-1}_{j_{dim}, i_{dim}} |J| \f]
   * where is jacobian of transformation from reference to physical coordinates at `i_qpoint`th quadrature point.
   * For Cartesian elements, returns `nullptr`.
   *
   * layout: [i_dim][j_dim][i_qpoint]
   */
  virtual double* reference_level_normals() = 0;
  //! \brief Jacobian determinant of reference -> physical coordinate transform
  //! \details for Cartesian elements, returns `nullptr`
  virtual double* jacobian_determinant() = 0;
  //! \brief where the extrapolated face area-weighted normal vectors are stored
  //! \details for Cartesian elements, returns `nullptr`
  virtual double* kernel_face_normal(int i_face) = 0;
};

}
#endif
