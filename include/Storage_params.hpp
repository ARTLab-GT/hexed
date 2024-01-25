#ifndef HEXED_STORAGE_PARAMS_HPP_
#define HEXED_STORAGE_PARAMS_HPP_

namespace hexed
{

//! \brief the parameters of the numerical scheme that are necessary to determine storage requirements
class Storage_params
{
  public:
  static constexpr int n_advection(int row_size) {return row_size;}

  int n_stage; //!< \brief number of time integration stages
  int n_var; //!< \brief number of independent physical state variables
  int n_dim; //!< \brief number of dimensions
  int row_size; //!< \brief \ref basis_row_size "row size" of basis
  int n_forcing = 4; //!< \brief number of artificial viscosity forcing variables

  int n_qpoint() const; //!< \brief number of quadrature points per element
  int n_dof() const; //!< \brief total number of physical degrees of freedom per element
  int size() const; //!< \brief total number of values that must be stored per element (including stages)
  int n_vertices() const; //!< \brief number of vertices of each element
  int n_var_numeric() const; //!< \brief number of _numerical_ variables (as opposed to physical state variables)
  int n_dof_numeric() const; //!< \brief number of numerical degrees of freedom
};

}
#endif
