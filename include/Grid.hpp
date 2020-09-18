#ifndef GRID_HPP_
#define GRID_HPP_

#include <vector>
#include <string>

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
  std::vector<double> state_r;
  std::vector<double> state_w;
  std::vector<int> pos;
  double mesh_size;
  double time;

  Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);
  virtual ~Grid();

  std::vector<double> get_pos(int i_elem);
  void visualize(std::string file_name);
  void print();

  private:
  void populate_slice(std::vector<double>&, std::vector<int>, int);
};

#endif
