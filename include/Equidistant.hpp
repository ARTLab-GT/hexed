#ifndef HEXED_EQUIDISTANT_HPP_
#define HEXED_EQUIDISTANT_HPP_

#include "Basis.hpp"

namespace hexed
{

/*! \brief `Basis` with equidistant quadrature points.
 * \details Used only for testing some basic featurs.
 * This would not be a good basis to use in practice
 * because of [Runge's phenomenon](https://en.wikipedia.org/wiki/Runge%27s_phenomenon).
 * As a result, some of the member functions not used in testing are placeholders we have not bothered to implement.
 */
class Equidistant : public Basis
{
  public:
  //! constructor
  Equidistant(int row_size_arg);
  double node(int i) const override;
  Eigen::MatrixXd diff_mat() const override;
  Eigen::VectorXd node_weights() const override; //!< placeholder (throws an exception)
  Eigen::MatrixXd boundary() const override;
  Eigen::VectorXd orthogonal(int degree) const override; //!< placeholder (throws an exception)
  Eigen::MatrixXd prolong(int i_half) const override; //!< placeholder (throws an exception)
  Eigen::MatrixXd restrict(int i_half) const override; //!< placeholder (throws an exception)
  double max_cfl_convective() const override; //!< placeholder (throws an exception)
  double max_cfl_diffusive() const override; //!< placeholder (throws an exception)
  double cancellation_convective() const override; //!< placeholder (throws an exception)
  double cancellation_diffusive() const override; //!< placeholder (throws an exception)
};

}
#endif
