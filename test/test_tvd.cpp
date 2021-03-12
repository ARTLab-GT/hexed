#include <catch.hpp>

#include <cartdgConfig.hpp>
#include <kernels/local/n_extrema.hpp>
#include <kernels/local/cpg_euler_matrix.hpp>
#include <kernels/local/cpg_euler_deformed.hpp>
#include <Gauss_lobatto.hpp>

TEST_CASE("n_extrema")
{
  double read [7] {};
  REQUIRE(cartdg::n_extrema<0>(read) == 0);
  REQUIRE(cartdg::n_extrema<1>(read) == 0);
  read[1] = 1.;
  REQUIRE(cartdg::n_extrema<2>(read) == 0);
  REQUIRE(cartdg::n_extrema<3>(read) == 1);
  read[1] = 0;
  REQUIRE(cartdg::n_extrema<3>(read) == 0);
  REQUIRE(cartdg::n_extrema<6>(read) == 0);
  read[1] = 1;
  read[2] = -1;
  REQUIRE(cartdg::n_extrema<6>(read) == 2);
  read[4] = -1;
  REQUIRE(cartdg::n_extrema<6>(read) == 4);
  read[3] = -2;
  REQUIRE(cartdg::n_extrema<6>(read) == 2);
}

TEST_CASE("regular TVD")
{
  constexpr int rank = MAX_BASIS_RANK;
  cartdg::Gauss_lobatto basis (rank);
  auto diff_mat = basis.diff_mat();
  auto weights = basis.node_weights();
  cartdg::Kernel_settings settings;
  settings.d_t_by_d_pos = 1.e-5;

  double read[4][rank][rank];
  double write[4][rank][rank] {};
  double mean [rank] {1., 1., 1.2, 2.5e5};
  for (int i_axis : {0, 1})
  {
    for (int i_var = 0; i_var < 4; ++i_var)
    {
      for (int i_row = 0; i_row < rank; ++i_row)
      {
        for (int j_row = 0; j_row < rank; ++j_row)
        {
          double mult = 1.;
          int row_inds [] {i_row, j_row};
          if ((row_inds[i_axis] > rank/2) && (i_var != 3))
          {
            mult -= 0.1;
          }
          read[i_var][i_row][j_row] = mean[i_var]*mult;
        }
      }
    }

    cartdg::cpg_euler_matrix<4, rank*rank, rank>(&read[0][0][0], &write[0][0][0], 1, diff_mat, weights, settings);
    double* var_read = &read[2][0][0];
    double* var_write = &write[2][0][0];
    for (int i_row = 0; i_row < rank; ++i_row)
    {
      double tv_read = 0.;
      double tv_write = 0.;
      for (int j_row = 1; j_row < rank; ++j_row)
      {
        int outer_stride = (i_axis == 0) ? 1 : rank;
        int inner_stride = (i_axis == 0) ? rank : 1;
        int ind0 = i_row*outer_stride + (j_row - 1)*inner_stride;
        int ind1 = i_row*outer_stride + j_row*inner_stride;
        tv_read += std::abs(var_read[ind1] - var_read[ind0]);
        tv_write += std::abs(var_write[ind1] - var_write[ind0]);
      }
      REQUIRE(tv_write <= tv_read);
    }
  }
}