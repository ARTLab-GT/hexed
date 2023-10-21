#ifndef HEXED_LINEAR_SYSTEM_HPP_
#define HEXED_LINEAR_SYSTEM_HPP_

#include "math.hpp"

namespace hexed
{

class Linear_equation
{
  public:
  virtual int n_vecs() = 0;
  virtual void scale(int output, int input, double scalar) = 0;
  virtual void add(int output, double coef0, int vec0, double coef1, int vec1) = 0;
  virtual double inner(int input0, int input1) = 0;
  double norm(int input);
  virtual void matvec(int output, int input) = 0;
};

class Dense_equation : public Linear_equation
{
  Mat<dyn, dyn> _mat;
  Mat<dyn, dyn> _vecs;
  public:
  Dense_equation(Mat<dyn, dyn> matrix, Mat<> rhs, int nv);
  Dense_equation(Mat<dyn, dyn> matrix, Mat<> rhs, int nv, Mat<> guess);
  inline auto vec(int i) {return _vecs(all, i);}
  int n_vecs() override;
  void scale(int output, int input, double scalar) override;
  void add(int output, double coef0, int vec0, double coef1, int vec1) override;
  double inner(int input0, int input1) override;
  void matvec(int output, int input) override;
};

}
#endif
