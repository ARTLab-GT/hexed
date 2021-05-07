#include <catch.hpp>

#include <cartdgConfig.hpp>
#include <get_local_cpg_euler.hpp>
#include <get_local_deformed_cpg_euler.hpp>
#include <get_req_visc_cpg_euler.hpp>
#include <get_av_flux.hpp>
#include <local/derivative.hpp>
#include <Gauss_lobatto.hpp>
#include <static_math.hpp>

class Identity_basis : public cartdg::Basis
{
  public:
  Identity_basis (int rank_arg) : Basis(rank_arg) {}
  double node(int i) { return 0.;}
  Eigen::MatrixXd diff_mat()
  {
    return -Eigen::MatrixXd::Identity(rank, rank);
  }
  Eigen::VectorXd node_weights()
  {
    return Eigen::VectorXd::Ones(rank);
  }
  Eigen::VectorXd orthogonal(int degree)
  {
    return Eigen::VectorXd::Zero(rank);
  }
};

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
    Identity_basis basis (rank);
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
    cartdg::get_local_cpg_euler(1, 2)(&read[0][0][0], &write[0][0][0], n_elem, basis, settings);
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
    Identity_basis basis (rank);
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
    cartdg::get_local_cpg_euler(3, 3)(&read[0][0][0], &write[0][0][0], n_elem, basis, settings);
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
    const int rank = std::min<int>(6, MAX_BASIS_RANK);
    double read [n_elem][4][rank*rank];
    double write[n_elem][4][rank*rank];
    cartdg::Gauss_lobatto basis (rank);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i = 0; i < rank; ++i)
      {
        for (int j = 0; j < rank; ++j)
        {
          int i_qpoint = i*rank + j;
          double mass = 1 + 0.1*basis.node(i) + 0.2*basis.node(j);
          double veloc0 = 10; double veloc1 = -20;
          double pres = 1e5*(1. - 0.3*basis.node(i) + 0.5*basis.node(j));
          double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1);

          read[i_elem][0][i_qpoint] = mass*veloc0;
          read[i_elem][1][i_qpoint] = mass*veloc1;
          read[i_elem][2][i_qpoint] = mass;
          read[i_elem][3][i_qpoint] = ener;
        }
      }
    }
    cartdg::get_local_cpg_euler(2, rank)(&read[0][0][0], &write[0][0][0], n_elem, basis, settings);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank*rank; ++i_qpoint)
      {
        REQUIRE((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint])
                == Approx(-0.1*(0.1*10 - 0.2*20)));
        REQUIRE((write[i_elem][3][i_qpoint] - read[i_elem][3][i_qpoint])
                == Approx(-0.1*(1e5*(1./0.4 + 1.)*(-0.3*10 - 0.5*20) + 0.5*(10*10 + 20*20)*(0.1*10 - 0.2*20))));
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
    Identity_basis basis (rank);
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
    cartdg::get_local_deformed_cpg_euler(1, 2)(&read[0][0][0], &write[0][0][0], 
                                               &jacobian[0][0][0][0], n_elem,
                                               basis, settings);
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
    Identity_basis basis (rank);
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
    cartdg::get_local_deformed_cpg_euler(3, 3)(&read[0][0][0], &write[0][0][0], 
                                               &jacobian[0][0][0][0], n_elem,
                                               basis, settings);
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
          double pres = 1e5*(1. - 0.3*pos0 + 0.5*pos1);
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
    cartdg::get_local_deformed_cpg_euler(2, rank)(&read[0][0][0], &write[0][0][0],
                                                  &jacobian[0][0][0][0], n_elem,
                                                  basis, settings);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < rank*rank; ++i_qpoint)
      {
        REQUIRE((write[i_elem][2][i_qpoint] - read[i_elem][2][i_qpoint])
                 == Approx(-0.1*(0.1*10 - 0.2*20)));
        REQUIRE((write[i_elem][3][i_qpoint] - read[i_elem][3][i_qpoint])
                == Approx(-0.1*(1e5*(1./0.4 + 1.)*(-0.3*10 - 0.5*20) + 0.5*(10*10 + 20*20)*(0.1*10 - 0.2*20))));
      }
    }
  }
}

TEST_CASE("derivative")
{
  const int row_size = MAX_BASIS_RANK;
  cartdg::Gauss_lobatto basis (row_size);
  cartdg::Kernel_settings settings;
  SECTION("calculations")
  {
    double read [row_size];
    double write [row_size];
    SECTION("constant function")
    {
      for (int i = 0; i < row_size; ++i)
      {
        read[i] = 1.;
        write[i] = 1.;
      }
      derivative<1, 1, row_size, row_size>(read, write, 1, 0, 0, 0, basis, settings);
      for (int i = 0; i < row_size; ++i)
      {
        REQUIRE(write[i] == Approx(0.).margin(1e-14));
      }
    }
    SECTION("polynomial")
    {
      for (int i = 0; i < row_size; ++i)
      {
        read[i] = std::pow(basis.node(i), 3);
        write[i] = 1.;
      }
      derivative<1, 1, row_size, row_size>(read, write, 1, 0, 0, 0, basis, settings);
      for (int i = 0; i < row_size; ++i)
      {
        REQUIRE(write[i] == Approx(3*std::pow(basis.node(i), 2)).margin(1e-14));
      }
    }
    SECTION("exponential") // Not mathematically exact but numerically very good
    {
      for (int i = 0; i < row_size; ++i)
      {
        read[i] = std::exp(basis.node(i));
        write[i] = 1.;
      }
      derivative<1, 1, row_size, row_size>(read, write, 1, 0, 0, 0, basis, settings);
      for (int i = 0; i < row_size; ++i)
      {
        REQUIRE(write[i] == Approx(std::exp(basis.node(i))).margin(1e-14));
      }
    }
  }
  SECTION("data juggling")
  {
    SECTION("3D")
    {
      const int n_qpoint = row_size*row_size*row_size;
      double coefs [] {1.103, -4.044, 0.392};
      for (int i_axis = 0; i_axis < 3; ++i_axis)
      {
        double read [row_size][row_size][row_size] {};
        double write [row_size][row_size][row_size] {};
        for (int i = 0; i < row_size; ++i)
        {
          for (int j = 0; j < row_size; ++j)
          {
            for (int k = 0; k < row_size; ++k)
            {
              int inds [] {i, j, k};
              for (int j_axis : {0, 1, 2}) read[i][j][k] += coefs[j_axis]*basis.node(inds[j_axis]);
            }
          }
        }
        derivative<1, 1, n_qpoint, row_size>(read[0][0], write[0][0], 1, 0, 0, i_axis, basis, settings);
        for (int i = 0; i < row_size; ++i)
        {
          for (int j = 0; j < row_size; ++j)
          {
            for (int k = 0; k < row_size; ++k)
            {
              REQUIRE(write[i][j][k] == Approx(coefs[i_axis]));
            }
          }
        }
      }
    }
    SECTION("multi-element")
    {
      double read [3][row_size] {};
      double write [3][row_size] {};
      double coefs [] {1.103, -4.044, 0.392};
      for (int i_elem = 0; i_elem < 3; ++i_elem)
      {
        for (int i = 0; i < row_size; ++i)
        {
          read[i_elem][i] = coefs[i_elem]*basis.node(i);
        }
      }
      derivative<1, 1, row_size, row_size>(read[0], write[0], 3, 0, 0, 0, basis, settings);
      for (int i_elem = 0; i_elem < 3; ++i_elem)
      {
        for (int i = 0; i < row_size; ++i)
        {
          REQUIRE(write[i_elem][i] == Approx(coefs[i_elem]));
        }
      }
    }
    SECTION("multivariable")
    {
      SECTION("1 var to 4 var")
      {
        double read [2][4][row_size] {};
        double write [2][row_size] {};
        for (int i_elem : {0, 1})
        {
          for (int i = 0; i < row_size; ++i)
          {
            read[i_elem][1][i] = std::pow(basis.node(i), 3);
          }
        }
        derivative<4, 1, row_size, row_size>(read[0][0], write[0], 2, 0, 0, 0, basis, settings);
        for (int i_elem : {0, 1})
        {
          for (int i = 0; i < row_size; ++i)
          {
            REQUIRE(write[i_elem][i] == Approx(0.).margin(1e-14));
          }
        }
        derivative<4, 1, row_size, row_size>(read[0][0], write[0], 2, 1, 0, 0, basis, settings);
        for (int i_elem : {0, 1})
        {
          for (int i = 0; i < row_size; ++i)
          {
            REQUIRE(write[i_elem][i] == Approx(3*std::pow(basis.node(i), 2)).margin(1e-14));
          }
        }
      }
      SECTION("3 var to 1 var")
      {
        double read [2][row_size] {};
        double write [2][3][row_size] {};
        for (int i_elem : {0, 1})
        {
          for (int i = 0; i < row_size; ++i)
          {
            read[i_elem][i] = std::pow(basis.node(i), 3);
          }
        }
        derivative<1, 3, row_size, row_size>(read[0], write[0][0], 2, 0, 1, 0, basis, settings);
        for (int i_elem : {0, 1})
        {
          for (int i = 0; i < row_size; ++i)
          {
            REQUIRE(write[i_elem][0][i] == Approx(0.).margin(1e-14));
          }
        }
        for (int i_elem : {0, 1})
        {
          for (int i = 0; i < row_size; ++i)
          {
            REQUIRE(write[i_elem][1][i] == Approx(3*std::pow(basis.node(i), 2)).margin(1e-14));
          }
        }
      }
    }
  }
}

TEST_CASE("req_visc")
{
  const int rank = MAX_BASIS_RANK;
  const int n_qpoint = cartdg::static_math::pow(rank, 3);
  cartdg::Gauss_lobatto basis (rank);
  cartdg::Kernel_settings settings;
  settings.d_pos = 0.5;
  double read [2][5][n_qpoint];
  double visc [2][2][2][2] {};
  for (int i_elem = 0; i_elem < 2; ++i_elem)
  {
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      read[i_elem][0][i_qpoint] = 0.;
      read[i_elem][1][i_qpoint] = 0.;
      read[i_elem][2][i_qpoint] = 0.;
      read[i_elem][3][i_qpoint] = 1.225;
      read[i_elem][4][i_qpoint] = 101325/0.4;
    }
  }
  read[1][3][3] = 1.5; // set an anomaly to trip the indicator in element 1
  cartdg::get_req_visc_cpg_euler(3, rank)(read[0][0], visc[0][0][0], 2, basis, settings);
  REQUIRE(visc[0][0][0][0] == 0.);
  REQUIRE(visc[0][1][1][1] == 0.);
  REQUIRE(visc[1][0][0][0] == Approx(0.5*340.29/(rank - 1.)).margin(0.01));
  REQUIRE(visc[1][1][1][1] == Approx(0.5*340.29/(rank - 1.)).margin(0.01));
}

TEST_CASE("av_flux")
{
  #if MAX_BASIS_RANK >= 3
  cartdg::Gauss_lobatto basis (3);
  cartdg::Kernel_settings settings;
  settings.d_pos = 0.5;
  settings.d_t_by_d_pos = 3.;
  double flux [2][9] {};
  double visc [2][4] {};
  visc[0][0] = 0.2;
  flux[1][1] = 0.3;
  flux[1][8] = 0.4;
  for (int i = 0; i < 4; ++i) visc[1][i] = 1.;
  for (int i = 0; i < 9; ++i) flux[0][i] = 1.;
  cartdg::get_av_flux(2, 3)(flux[0], visc[0], 2, basis, settings);

  REQUIRE(flux[0][0] == Approx(0.2*1.5));
  REQUIRE(flux[0][1] == Approx(0.1*1.5));
  REQUIRE(flux[0][2] == Approx(0.0*1.5));
  REQUIRE(flux[0][3] == Approx(0.1*1.5));
  REQUIRE(flux[0][4] == Approx(0.05*1.5));
  REQUIRE(flux[0][5] == Approx(0.0*1.5));

  REQUIRE(flux[1][0] == Approx(0.0));
  REQUIRE(flux[1][1] == Approx(0.3*1.5));
  REQUIRE(flux[1][8] == Approx(0.4*1.5));
  #endif
}
