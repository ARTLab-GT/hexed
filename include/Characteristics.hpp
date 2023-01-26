#ifndef HEXED_CHARACTERISTICS_HPP_
#define HEXED_CHARACTERISTICS_HPP_

#include "math.hpp"

namespace hexed
{

class Characteristics
{
  Mat<3> vals;
  Mat<3, 3> vecs;
  Mat<> dir;
  const int n_dim;
  double mass;
  Mat<> veloc;
  double nrml(Mat<> vec);
  Mat<> tang(Mat<> vec);
  public:
  Characteristics(Mat<> state, Mat<> direction);
  inline Mat<3> eigvals() {return vals;}
  Mat<dyn, 3> decomp(Mat<> state);
};

}
#endif
