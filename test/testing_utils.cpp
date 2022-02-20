#include <cartdgConfig.hpp>
#include "testing_utils.hpp"

std::unique_ptr<cartdg::Solution> deformed_test_setup_2d(bool perturb_face)
{
  // I know this setup is complicated... sorry :(
  const int row_size = cartdg::config::max_row_size;
  std::unique_ptr<cartdg::Solution> sol_ptr {new cartdg::Solution {4, 2, row_size, 1.}};
  cartdg::Solution& sol = *sol_ptr;
  sol.add_empty_grid(1);
  cartdg::Regular_grid& grid = sol.reg_grids[0];
  sol.add_deformed_grid(1);
  cartdg::Deformed_grid& def_grid = sol.def_grids[0];

  /*
   * Element diagram.
   * [i] = grid.element(i)
   * /i/ = def_grid.deformed_element(i)
   *
   * x1  ^
   * 2.5 _| [ 5][16][17][18][19][11]
   * 2.0 _| [ 4]/ 3// 7//11//15/[10]
   * 1.5 _| [ 3]/ 2// 6//10//14/[ 9]
   * 1.0 _| [ 2]/ 1// 5// 9//13/[ 8]
   * 0.5 _| [ 1]/ 0// 4// 8//12/[ 7]
   * 0.0 _| [ 0][12][13][14][15][ 6]
   *      + ------------------------>
   *        |   |   |   |   |   |
   *        0.0 0.5 1.0 1.5 2.0 2.5 x0
   */
  for (int i = 0; i < 6; ++i) grid.add_element({0, i});
  for (int i = 0; i < 6; ++i) grid.add_element({5, i});
  for (int i = 0; i < 4; ++i) grid.add_element({i + 1, 0});
  for (int i = 0; i < 4; ++i) grid.add_element({i + 1, 5});
  grid.auto_connect({6, 6});

  for (int i = 1; i < 5; ++i) {
    for (int j = 1; j < 5; ++j) {
      def_grid.add_element({i, j});
    }
  }
  double center [] {1.5 + .1, 1.5 + .1};
  {
    // move one vertex to create deformation
    cartdg::Deformed_element& elem {def_grid.deformed_element(5)};
    elem.vertex(3).pos = {center[0], center[1], 0.};
    if (perturb_face) elem.node_adjustments()[row_size + 1] = 0.1;
  }
  // deform corresponding vertices to match
  {
    cartdg::Deformed_element& elem {def_grid.deformed_element(9)};
    elem.vertex(1).pos = {center[0], center[1], 0.};
  }
  {
    cartdg::Deformed_element& elem {def_grid.deformed_element(10)};
    elem.vertex(0).pos = {center[0], center[1], 0.};
  }
  {
    // rotate element to facilitate testing connections along different axes
    cartdg::Deformed_element& elem {def_grid.deformed_element(6)};
    elem.vertex(0).pos = {1.0, 2.0, 0.};
    elem.vertex(1).pos = {1.5, 2.0, 0.};
    elem.vertex(2).pos = {1.0, 1.5, 0.};
    elem.vertex(3).pos = {center[0], center[1], 0.};
  }
  {
    // perturb another vertex (don't bother to move corresponding vertices)
    cartdg::Deformed_element& elem {def_grid.deformed_element(15)};
    elem.vertex(0).pos[0] += 0.2;
    elem.vertex(0).pos[1] += 0.2;
  }
  // connect all the deformed elements except 6 (the rotated one)
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 4; ++j) {
      int inds [2][2] {{4*i+j, 4*(i+1)+j}, {4*j+i, 4*j+i+1}};
      for (int i_dim : {0, 1}) {
        if ((inds[i_dim][0] != 6) && (inds[i_dim][1] != 6)) {
          def_grid.connect({inds[i_dim][0], inds[i_dim][1]}, {i_dim, i_dim}, {1, 0});
        }
      }
    }
  }
  // connect element 6
  def_grid.connect({6, 7}, {0, 1}, {0, 0});
  def_grid.connect({6, 5}, {0, 1}, {1, 1});
  def_grid.connect({6,  2}, {1, 0}, {0, 1});
  def_grid.connect({6, 10}, {1, 0}, {1, 0});
  // connect deformed to non-deformed
  for (int i = 0; i < 4; ++i) {
    grid.add_connection(&grid.element(i + 1), &def_grid.element(i), 0);
    grid.add_connection(&def_grid.element(i + 12), &grid.element(i + 7), 0);
    grid.add_connection(&grid.element(i + 12), &def_grid.element(4*i), 1);
    grid.add_connection(&def_grid.element(4*i + 3), &grid.element(i + 16), 1);
  }
  def_grid.calc_jacobian();
  return sol_ptr;
}
