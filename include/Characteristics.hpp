#ifndef HEXED_CHARACTERISTICS_HPP_
#define HEXED_CHARACTERISTICS_HPP_

#include "math.hpp"

namespace hexed
{

/*
 * Decomposes state vectors into characteristics
 * which are eigenvectors of the Jacobian of the inviscid flux function.
 * This is useful for characteristic-based boundary conditions.
 */
class Characteristics
{
  // 1D eigensystem
  Mat<3> vals;
  Mat<3, 3> vecs;
  Mat<> dir; // normalized flux direction
  const int n_dim;
  // some properties of the reference state
  double mass;
  Mat<> veloc;
  double nrml(Mat<> vec); // compute component of vector normal to `dir`
  Mat<> tang(Mat<> vec); // remove normal component of vector to obtain part tangential to `dir`

  public:
  // construct with a direction in which to compute the flux
  // and a reference state vector about which to compute the Jacobian
  Characteristics(Mat<> state, Mat<> direction);
  // get eigenvalues of Jacobian
  inline Mat<3> eigvals() {return vals;}
  /*
   * Decompose a state vector into eigenspaces.
   * Column `j` should be an eigenvector of the Jacobian with eigenvalue `eigvals()(j)`
   * and the sum of the columns should be `state`.
   */
  Mat<dyn, 3> decomp(Mat<> state);
};

}
#endif
