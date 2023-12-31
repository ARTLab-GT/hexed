#ifndef HEXED_BASIS_HPP_
#define HEXED_BASIS_HPP_

#include <Eigen/Dense>

namespace hexed
{

/*! \brief Represents a basis of [Lagrange polynomials](https://en.wikipedia.org/wiki/Lagrange_polynomial).
 * \details Interpolates between a set of quadrature points, a.k.a. nodes.
 * The `i`th basis vector is the polynomial with a value
 * of 1 at the `i`th node and 0 at all other nodes.
 */
class Basis
{
  protected:
  //! \brief nondimensional minimum real part of eigenvals of 1D convection operator
  virtual double min_eig_convection() const = 0;
  virtual double quadratic_safety() const = 0;

  public:
  const int row_size; //!< \brief \ref basis_row_size "row size"

  //! \param row_size_arg specifies member `row_size`.
  Basis(int row_size_arg);
  virtual ~Basis();
  //! \brief the `i`th interpolation node (i.e. quadrature point). Should be in interval [0, 1].
  virtual double node(int i) const = 0;
  Eigen::VectorXd nodes() const; //!< \brief Get vector of nodes
  //! \brief `node_weights()(i)` is the quadrature weight associated with `node(i)` Integrates in domain [0, 1].
  virtual Eigen::VectorXd node_weights() const = 0;
  //! \brief Differentiation matrix. `diff_mat()(i, j)` is the derivative of interpolating polynomial `j` evaluated at `node(i)`.
  virtual Eigen::MatrixXd diff_mat() const = 0;
  //! \brief `boundary()(i, j)` is the `j`th basis polynomial evaluated at `i` for `i` in {0, 1}.
  virtual Eigen::MatrixXd boundary() const = 0;
  /*! \brief `orthogonal(degree)(i)` is the `degree`th-degree [Legendre polynomial](https://en.wikipedia.org/wiki/Legendre_polynomials) evaluated at `node(i)`.
   * \details Polynomial is scaled so that its norm is 1 with respect to the quadrature rule of the basis.
   */
  virtual Eigen::VectorXd orthogonal(int degree) const = 0;
  //! \brief Matrix that applies a low-pass filter to the Legendre modes of a polynomial
  virtual Eigen::MatrixXd filter() const = 0;
  /*! \brief Interpolate to arbitrary points.
   * \details Current implementation is not very performance-optimized. `interpolate(sample)(i, j)` is the `j`th basis polynomial evaluated at `sample(i)`
   */
  Eigen::MatrixXd interpolate(const Eigen::VectorXd& sample) const;
  /*! \brief Matrix to prolong into polynomial space of one higher refinement level.
   * \details `prolong(i_half)(i, j)` is the `j`th basis polynomial evaluated at the `i`th node in the refined space.
   * If `i_half` is 0, refined space is interval [0, 0.5]. if `i_half` is 1, then [0.5, 1]
   */
  virtual Eigen::MatrixXd prolong(int i_half) const = 0;
  /*! \brief Restrict from refined space above to polynomial space spanned by this basis.
   * \details `restrict(i_half)(i, j)` is the value at the `i`th standard node
   * of the orthogonal projection of the `j`th polynomial in the
   * refined space into the space of this basis.
   */
  virtual Eigen::MatrixXd restrict(int i_half) const = 0;
  //! \brief maximum stable CFL number for 1D convection
  double max_cfl() const;
  //! \brief ratio of time step for convection stage 2 to stage 1
  double step_ratio() const;
  //! \brief nondimensional minimum eigenvalue of 1D diffusion operator
  virtual double min_eig_diffusion() const = 0;
};

}
#endif
