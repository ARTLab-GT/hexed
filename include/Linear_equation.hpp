#ifndef HEXED_LINEAR_SYSTEM_HPP_
#define HEXED_LINEAR_SYSTEM_HPP_

#include "math.hpp"

namespace hexed
{

/*! \brief represents a linear equation in an arbitrary vector space with a black-box operator
 * \details represents an equation of the form \f$ Ax = b \f$.
 * This object manages (or has access to) storage for some number of vectors indexed 0 to `n_vecs()`.
 * The member functions of this class provide basic vector operations including addition, scalar multiplication, and inner products
 * as well as `matvec()` to compute matrix-vector products.
 * Derived classes may also implement other ways to access the vector values.
 * Vector 0 is the current best estimate for the solution \f$ x \f$ and vector 1 is the right hand side \f$ b \f$.
 * All other vectors are available for solver algorithms to use.
 * \note In the member functions below, any number of the output and/or input vectors may be the same.
 * Derived classes must be implemented in such a way that this does not create aliasing problems.
 * \note This clas was created in order to experiment with implicit solvers.
 * It is not used in any of the production-level code to date.
 */
class Linear_equation
{
  public:
  virtual int n_vecs() = 0; //!< \brief number of vectors this object provides access to (including 0 and 1)
  virtual void scale(int output, int input, double scalar) = 0; //!< \brief scales vector `input` by `scalar` and writes the result to vector `output`
  //! \brief assigns \f$ z = ax + by \f$ where \f$ z, a, x, b, y \f$ are given by the arguments in order
  virtual void add(int output, double coef0, int vec0, double coef1, int vec1) = 0;
  virtual double inner(int input0, int input1) = 0; //!< \brief returns the inner product of the two vectors specified by the inputs
  double norm(int input); //!< \brief norm induced by `inner()`
  virtual void matvec(int output, int input) = 0; //!< \brief computes the product of the operator \f$ A \f$ with vector `input` and writes it to `output`.
};

//! \brief equation in \f$ \mathbb{R}^n \f$ with a dense matrix
class Dense_equation : public Linear_equation
{
  Mat<dyn, dyn> _mat;
  Mat<dyn, dyn> _vecs;
  public:
  //! \brief constructs an equation from \f$ A \f$ and \f$ b \f$ with an optional initial guess (otherwise defaults to the zero vector)
  Dense_equation(Mat<dyn, dyn> matrix, Mat<> rhs, int nv);
  Dense_equation(Mat<dyn, dyn> matrix, Mat<> rhs, int nv, Mat<> guess); //!< \overload
  inline auto vec(int i) {return _vecs(all, i);} //!< gets a view of a vector as an Eigen object
  int n_vecs() override;
  void scale(int output, int input, double scalar) override;
  void add(int output, double coef0, int vec0, double coef1, int vec1) override;
  double inner(int input0, int input1) override;
  void matvec(int output, int input) override;
};

}
#endif
