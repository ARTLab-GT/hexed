#include <catch.hpp>
#include <kernels/local/cpg_euler_matrix.hpp>
#include <kernels/local/cpg_euler_deformed.hpp>

#include <Gauss_lobatto.hpp>

TEST_CASE("CPG Euler matrix form")
{
  cartdg::Kernel_settings settings;
  settings.d_t_by_d_pos = 0.1;
  SECTION("1D")
  {
    const int n_elem = 5;
    const int rank = 2;
    double read [n_elem][3][rank];
    double write[n_elem][3][rank];
    Eigen::MatrixXd diff_mat = -Eigen::MatrixXd::Identity(rank, rank);
    double mass = 1.225; double veloc = 10; double pres = 1e5;
    double mmtm = mass*veloc; double ener = pres/0.4 + 0.5*mass*veloc*veloc;
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank; ++i_qpoint)
      {
          read[i_elem][0][i_qpoint] = mmtm;
          read[i_elem][1][i_qpoint] = mass;
          read[i_elem][2][i_qpoint] = ener;
      }
    }
    cartdg::cpg_euler_matrix<3, 2, 2>(&read[0][0][0], &write[0][0][0], n_elem,
                                      diff_mat, settings);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank; ++i_qpoint)
      {
        REQUIRE((write[i_elem][0][i_qpoint] - read[i_elem][0][i_qpoint])
              == Approx(0.1*(mmtm*mmtm/mass + pres)));
        REQUIRE((write[i_elem][1][i_qpoint] - read[i_elem][1][i_qpoint])
              == Approx(0.1*mmtm));
        REQUIRE((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint])
              == Approx(0.1*((ener + pres)*veloc)));
      }
    }
  }

  SECTION("3D")
  {
    const int n_elem = 5;
    const int rank = 3;
    double read [n_elem][5][rank*rank*rank];
    double write[n_elem][5][rank*rank*rank];
    Eigen::MatrixXd diff_mat = -Eigen::MatrixXd::Identity(rank, rank);
    double mass = 1.225;
    double veloc0 = 10; double veloc1 = -20; double veloc2 = 30;
    double pres = 1e5;
    double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1 + veloc2*veloc2);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank*rank*rank; ++i_qpoint)
      {
          read[i_elem][0][i_qpoint] = mass*veloc0;
          read[i_elem][1][i_qpoint] = mass*veloc1;
          read[i_elem][2][i_qpoint] = mass*veloc2;
          read[i_elem][3][i_qpoint] = mass;
          read[i_elem][4][i_qpoint] = ener;
      }
    }
    cartdg::cpg_euler_matrix<5, 27, 3>(&read[0][0][0], &write[0][0][0], n_elem,
                                       diff_mat, settings);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank*rank*rank; ++i_qpoint)
      {
        REQUIRE((write[i_elem][0][i_qpoint] - read[i_elem][0][i_qpoint])
              == Approx(0.1*(mass*veloc0*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE((write[i_elem][1][i_qpoint] - read[i_elem][1][i_qpoint])
              == Approx(0.1*(mass*veloc1*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint])
              == Approx(0.1*(mass*veloc2*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE((write[i_elem][3][i_qpoint] - read[i_elem][3][i_qpoint])
              == Approx(0.1*mass*(veloc0 + veloc1 + veloc2)));
        REQUIRE((write[i_elem][4][i_qpoint] - read[i_elem][4][i_qpoint])
              == Approx(0.1*((ener + pres)*(veloc0 + veloc1 + veloc2))));
      }
    }
  }

  SECTION("2D non-constant")
  {
    const int n_elem = 5;
    const int rank = 6;
    double read [n_elem][4][rank*rank];
    double write[n_elem][4][rank*rank];
    cartdg::Gauss_lobatto basis (rank);
    Eigen::MatrixXd diff_mat = basis.diff_mat();
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i = 0; i < rank; ++i)
      {
        for (int j = 0; j < rank; ++j)
        {
          int i_qpoint = i*rank + j;
          double mass = 1 + 0.1*basis.node(i) + 0.2*basis.node(j);
          double veloc0 = 10; double veloc1 = -20;
          double pres = 1e5;
          double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1);

          read[i_elem][0][i_qpoint] = mass*veloc0;
          read[i_elem][1][i_qpoint] = mass*veloc1;
          read[i_elem][2][i_qpoint] = mass;
          read[i_elem][3][i_qpoint] = ener;
        }
      }
    }
    cartdg::cpg_euler_matrix<4, rank*rank, rank>(&read[0][0][0], &write[0][0][0], n_elem,
                                                 diff_mat, settings);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank*rank; ++i_qpoint)
      {
        REQUIRE((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint])
                == Approx(-0.1*(0.1*10 - 0.2*20)));
      }
    }
  }
}

TEST_CASE("CPG Euler deformed elements")
{
  cartdg::Kernel_settings settings;
  settings.d_t_by_d_pos = 0.1;
  SECTION("1D regular")
  {
    const int n_elem = 5;
    const int rank = 2;
    double read [n_elem][3][rank];
    double write[n_elem][3][rank];
    double jacobian [n_elem][1][1][rank] {};
    Eigen::MatrixXd diff_mat = -Eigen::MatrixXd::Identity(rank, rank);
    double mass = 1.225; double veloc = 10; double pres = 1e5;
    double mmtm = mass*veloc; double ener = pres/0.4 + 0.5*mass*veloc*veloc;
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank; ++i_qpoint)
      {
          read[i_elem][0][i_qpoint] = mmtm;
          read[i_elem][1][i_qpoint] = mass;
          read[i_elem][2][i_qpoint] = ener;
          jacobian[i_elem][0][0][i_qpoint] = 1.;
      }
    }
    cartdg::cpg_euler_deformed<3, 2, 2>(&read[0][0][0], &write[0][0][0], 
                                        &jacobian[0][0][0][0], n_elem,
                                        diff_mat, settings);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank; ++i_qpoint)
      {
        REQUIRE((write[i_elem][0][i_qpoint] - read[i_elem][0][i_qpoint])
              == Approx(0.1*(mmtm*mmtm/mass + pres)));
        REQUIRE((write[i_elem][1][i_qpoint] - read[i_elem][1][i_qpoint])
              == Approx(0.1*mmtm));
        REQUIRE((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint])
              == Approx(0.1*((ener + pres)*veloc)));
      }
    }
  }

  SECTION("3D regular")
  {
    const int n_elem = 5;
    const int rank = 3;
    double read [n_elem][5][rank*rank*rank];
    double write[n_elem][5][rank*rank*rank];
    double jacobian [n_elem][3][3][rank*rank*rank] {};
    Eigen::MatrixXd diff_mat = -Eigen::MatrixXd::Identity(rank, rank);
    double mass = 1.225;
    double veloc0 = 10; double veloc1 = -20; double veloc2 = 30;
    double pres = 1e5;
    double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1 + veloc2*veloc2);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank*rank*rank; ++i_qpoint)
      {
          read[i_elem][0][i_qpoint] = mass*veloc0;
          read[i_elem][1][i_qpoint] = mass*veloc1;
          read[i_elem][2][i_qpoint] = mass*veloc2;
          read[i_elem][3][i_qpoint] = mass;
          read[i_elem][4][i_qpoint] = ener;
          jacobian[i_elem][0][0][i_qpoint] = 1.;
          jacobian[i_elem][1][1][i_qpoint] = 1.;
          jacobian[i_elem][2][2][i_qpoint] = 1.;
      }
    }
    cartdg::cpg_euler_deformed<5, 27, 3>(&read[0][0][0], &write[0][0][0], 
                                         &jacobian[0][0][0][0], n_elem,
                                         diff_mat, settings);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank*rank*rank; ++i_qpoint)
      {
        REQUIRE((write[i_elem][0][i_qpoint] - read[i_elem][0][i_qpoint])
              == Approx(0.1*(mass*veloc0*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE((write[i_elem][1][i_qpoint] - read[i_elem][1][i_qpoint])
              == Approx(0.1*(mass*veloc1*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint])
              == Approx(0.1*(mass*veloc2*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE((write[i_elem][3][i_qpoint] - read[i_elem][3][i_qpoint])
              == Approx(0.1*mass*(veloc0 + veloc1 + veloc2)));
        REQUIRE((write[i_elem][4][i_qpoint] - read[i_elem][4][i_qpoint])
              == Approx(0.1*((ener + pres)*(veloc0 + veloc1 + veloc2))));
      }
    }
  }

  SECTION("2D non-constant deformed")
  {
    const int n_elem = 5;
    const int rank = 6;
    double read [n_elem][4][rank*rank];
    double write[n_elem][4][rank*rank];
    double jacobian [n_elem][2][2][rank*rank] {};
    cartdg::Gauss_lobatto basis (rank);
    Eigen::MatrixXd diff_mat = basis.diff_mat();
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i = 0; i < rank; ++i)
      {
        for (int j = 0; j < rank; ++j)
        {
          int i_qpoint = i*rank + j;
          double pos0 = basis.node(i)*(1. - 0.5*basis.node(j));
          double pos1 = basis.node(j)*(1. - 0.3*basis.node(i));
          double mass = 1 + 0.1*pos0 + 0.2*pos1;
          double veloc0 = 10; double veloc1 = -20;
          double pres = 1e5;
          double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1);

          read[i_elem][0][i_qpoint] = mass*veloc0;
          read[i_elem][1][i_qpoint] = mass*veloc1;
          read[i_elem][2][i_qpoint] = mass;
          read[i_elem][3][i_qpoint] = ener;

          jacobian[i_elem][0][0][i_qpoint] = 1. - 0.5*basis.node(j);
          jacobian[i_elem][0][1][i_qpoint] = -0.5*basis.node(i);
          jacobian[i_elem][1][0][i_qpoint] = -0.3*basis.node(j);
          jacobian[i_elem][1][1][i_qpoint] = 1. - 0.3*basis.node(i);
        }
      }
    }
    cartdg::cpg_euler_deformed<4, rank*rank, rank>(&read[0][0][0], &write[0][0][0],
                                                   &jacobian[0][0][0][0], n_elem,
                                                   diff_mat, settings);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank*rank; ++i_qpoint)
      {
        CHECK((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint])
                == Approx(-0.1*(0.1*10 - 0.2*20)));
      }
    }
  }
}
