#ifndef GRID_HPP_
#define GRID_HPP_

#include <vector>

#include "Basis.hpp"

class Grid
{
  public:
  Basis& basis;
  int n_var;
  int n_dim;
  int n_qpoint;
  int n_dof;
  int n_elem;
  double * access_r;
  double * access_w;

  Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, Basis& basis_arg);
  virtual ~Grid();

  protected:
  std::vector<double> storage_r;
  std::vector<double> storage_w;
};

#endif
