#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <get_local_convective.hpp>
#include <get_local_deformed_convective.hpp>
#include <get_req_visc_regular_convective.hpp>
#include <get_av_flux.hpp>
#include <local/variable_derivative.hpp>
#include <get_local_derivative.hpp>
#include <Gauss_lobatto.hpp>
#include <math.hpp>

class Identity_basis : public cartdg::Basis
{
  public:
  Identity_basis (int row_size_arg) : Basis(row_size_arg) {}
  double node(int i) { return 0.;}
  Eigen::MatrixXd diff_mat()
  {
    return -Eigen::MatrixXd::Identity(row_size, row_size);
  }
  Eigen::VectorXd node_weights()
  {
    return Eigen::VectorXd::Ones(row_size);
  }
  Eigen::MatrixXd boundary()
  {
    return Eigen::MatrixXd::Zero(2, row_size);
  }
  Eigen::VectorXd orthogonal(int degree)
  {
    return Eigen::VectorXd::Zero(row_size);
  }
};

TEST_CASE("Local convective")
{
  cartdg::Kernel_settings settings;
  settings.d_t_by_d_pos = 0.1;
  SECTION("1D")
  {
    unsigned n_elem = 5;
    cartdg::Storage_params params {2, 3, 1, 2};
    unsigned n_qpoint = params.n_qpoint();
    cartdg::elem_vec elements;
    Identity_basis basis {int(params.row_size)};
    double mass = 1.225; double veloc = 10; double pres = 1e5;
    double mmtm = mass*veloc; double ener = pres/0.4 + 0.5*mass*veloc*veloc;
    for (unsigned i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new cartdg::Element {params});
      double* read = elements[i_elem]->stage(settings.i_read);
      for (unsigned i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
          read[0*n_qpoint + i_qpoint] = mmtm;
          read[1*n_qpoint + i_qpoint] = mass;
          read[2*n_qpoint + i_qpoint] = ener;
      }
    }
    cartdg::get_local_convective(1, 2)(elements, basis, settings);
    for (auto& element : elements)
    {
      double* r = element->stage(settings.i_read);
      double* w = element->stage(settings.i_write);
      for (unsigned i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        REQUIRE(w[0*n_qpoint + i_qpoint] - r[0*n_qpoint + i_qpoint] == Approx(0.1*(mmtm*mmtm/mass + pres)));
        REQUIRE(w[1*n_qpoint + i_qpoint] - r[1*n_qpoint + i_qpoint] == Approx(0.1*mmtm));
        REQUIRE(w[2*n_qpoint + i_qpoint] - r[2*n_qpoint + i_qpoint] == Approx(0.1*((ener + pres)*veloc)));
      }
    }
  }

  SECTION("3D")
  {
    const unsigned n_elem = 5;
    cartdg::Storage_params params {2, 5, 3, 3};
    unsigned n_qpoint = params.n_qpoint();
    cartdg::elem_vec elements;
    Identity_basis basis (params.row_size);
    double mass = 1.225;
    double veloc0 = 10; double veloc1 = -20; double veloc2 = 30;
    double pres = 1e5;
    double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1 + veloc2*veloc2);
    for (unsigned i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new cartdg::Element {params});
      double* read = elements[i_elem]->stage(settings.i_read);
      for (unsigned i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
          read[0*n_qpoint + i_qpoint] = mass*veloc0;
          read[1*n_qpoint + i_qpoint] = mass*veloc1;
          read[2*n_qpoint + i_qpoint] = mass*veloc2;
          read[3*n_qpoint + i_qpoint] = mass;
          read[4*n_qpoint + i_qpoint] = ener;
      }
    }
    cartdg::get_local_convective(3, 3)(elements, basis, settings);
    for (auto& element : elements)
    {
      double* r = element->stage(settings.i_read);
      double* w = element->stage(settings.i_write);
      for (unsigned i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        REQUIRE(w[0*n_qpoint + i_qpoint] - r[0*n_qpoint + i_qpoint]
                == Approx(0.1*(mass*veloc0*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE(w[1*n_qpoint + i_qpoint] - r[1*n_qpoint + i_qpoint]
                == Approx(0.1*(mass*veloc1*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE(w[2*n_qpoint + i_qpoint] - r[2*n_qpoint + i_qpoint]
                == Approx(0.1*(mass*veloc2*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE(w[3*n_qpoint + i_qpoint] - r[3*n_qpoint + i_qpoint]
                == Approx(0.1*mass*(veloc0 + veloc1 + veloc2)));
        REQUIRE(w[4*n_qpoint + i_qpoint] - r[4*n_qpoint + i_qpoint]
                == Approx(0.1*((ener + pres)*(veloc0 + veloc1 + veloc2))));
      }
    }
  }

  SECTION("2D non-constant")
  {
    const unsigned n_elem = 5;
    const unsigned row_size = std::min<unsigned>(6, unsigned(CARTDG_MAX_BASIS_ROW_SIZE));
    cartdg::Storage_params params {3, 4, 2, row_size};
    unsigned int n_qpoint = params.n_qpoint();
    cartdg::elem_vec elements;
    cartdg::Gauss_lobatto basis (row_size);
    for (unsigned i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new cartdg::Element {params});
      double* read = elements[i_elem]->stage(settings.i_read);
      for (unsigned i = 0; i < row_size; ++i)
      {
        for (unsigned j = 0; j < row_size; ++j)
        {
          int i_qpoint = i*row_size + j;
          double mass = 1 + 0.1*basis.node(i) + 0.2*basis.node(j);
          double veloc0 = 10; double veloc1 = -20;
          double pres = 1e5*(1. - 0.3*basis.node(i) + 0.5*basis.node(j));
          double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1);

          read[0*n_qpoint + i_qpoint] = mass*veloc0;
          read[1*n_qpoint + i_qpoint] = mass*veloc1;
          read[2*n_qpoint + i_qpoint] = mass;
          read[3*n_qpoint + i_qpoint] = ener;
        }
      }
    }
    cartdg::get_local_convective(2, row_size)(elements, basis, settings);
    for (auto& element : elements)
    {
      double* r = element->stage(settings.i_read);
      double* w = element->stage(settings.i_write);
      for (unsigned i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        REQUIRE(w[2*n_qpoint + i_qpoint] - r[2*n_qpoint + i_qpoint]
                == Approx(-0.1*(0.1*10 - 0.2*20)));
        REQUIRE(w[3*n_qpoint + i_qpoint] - r[3*n_qpoint + i_qpoint]
                == Approx(-0.1*(1e5*(1./0.4 + 1.)*(-0.3*10 - 0.5*20)
                          + 0.5*(10*10 + 20*20)*(0.1*10 - 0.2*20))));
      }
    }
  }
}

TEST_CASE("CPG Euler deformed elements")
{
  cartdg::Kernel_settings settings;
  settings.d_t_by_d_pos = 0.1;
  #if 0
  SECTION("1D regular")
  {
    const int n_elem = 5;
    const int row_size = 2;
    double read [n_elem][3][row_size];
    double write[n_elem][3][row_size];
    double jacobian [n_elem][1][1][row_size] {};
    Identity_basis basis (row_size);
    double mass = 1.225; double veloc = 10; double pres = 1e5;
    double mmtm = mass*veloc; double ener = pres/0.4 + 0.5*mass*veloc*veloc;
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
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
      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
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
    const int row_size = 3;
    double read [n_elem][5][row_size*row_size*row_size];
    double write[n_elem][5][row_size*row_size*row_size];
    double jacobian [n_elem][3][3][row_size*row_size*row_size] {};
    Identity_basis basis (row_size);
    double mass = 1.225;
    double veloc0 = 10; double veloc1 = -20; double veloc2 = 30;
    double pres = 1e5;
    double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1 + veloc2*veloc2);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < row_size*row_size*row_size; ++i_qpoint)
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
      for (int i_qpoint = 0; i_qpoint < row_size*row_size*row_size; ++i_qpoint)
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
  #endif

  SECTION("Local convective for deformed elements (2D non-constant)")
  {
    const int n_elem = 5;
    const int row_size = 6;
    cartdg::Gauss_lobatto basis (row_size);
    cartdg::Storage_params params {3, 4, 2, row_size};
    int n_qpoint = params.n_qpoint();
    cartdg::def_elem_vec elements;
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new cartdg::Deformed_element {params});
      double* stage    = elements.back()->stage(0);
      double* jacobian = elements.back()->jacobian();
      for (int i = 0; i < row_size; ++i)
      {
        for (int j = 0; j < row_size; ++j)
        {
          int i_qpoint = i*row_size + j;
          double pos0 = basis.node(i)*(1. - 0.5*basis.node(j));
          double pos1 = basis.node(j)*(1. - 0.3*basis.node(i));
          double mass = 1 + 0.1*std::pow(pos0, 5) - 0.3*std::pow(pos0, 2)*std::pow(pos1, 3) + 0.2*pos1;
          double veloc0 = 10; double veloc1 = -20;
          double pres = 1e5*(1. - 0.3*pos0 + 0.5*pos1);
          double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1);

          stage[0*n_qpoint + i_qpoint] = mass*veloc0;
          stage[1*n_qpoint + i_qpoint] = mass*veloc1;
          stage[2*n_qpoint + i_qpoint] = mass;
          stage[3*n_qpoint + i_qpoint] = ener;

          jacobian[(0*2 + 0)*n_qpoint + i_qpoint] = 1. - 0.5*basis.node(j);
          jacobian[(0*2 + 1)*n_qpoint + i_qpoint] = -0.5*basis.node(i);
          jacobian[(1*2 + 0)*n_qpoint + i_qpoint] = -0.3*basis.node(j);
          jacobian[(1*2 + 1)*n_qpoint + i_qpoint] = 1. - 0.3*basis.node(i);
        }
      }
    }

    cartdg::get_local_deformed_convective(2, row_size)(elements, basis, settings);

    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      double* stage_r  = elements[i_elem]->stage(0);
      double* stage_w  = elements[i_elem]->stage(1);
      for (int i = 0; i < row_size; ++i)
      {
        for (int j = 0; j < row_size; ++j)
        {
          int i_qpoint = i*row_size + j;
          double pos0 = basis.node(i)*(1. - 0.5*basis.node(j));
          double pos1 = basis.node(j)*(1. - 0.3*basis.node(i));
          double veloc_dot_grad_mass = (0.1*5*std::pow(pos0, 4)
                                        - 0.3*2*std::pow(pos0, 1)*std::pow(pos1, 3))*10
                                        + (-0.3*std::pow(pos0, 2)*3*std::pow(pos1, 2) + 0.2)*-20;
          REQUIRE((stage_w[2*n_qpoint + i_qpoint] - stage_r[2*n_qpoint + i_qpoint])
                   == Approx(-0.1*veloc_dot_grad_mass));
          REQUIRE((stage_w[3*n_qpoint + i_qpoint] - stage_r[3*n_qpoint + i_qpoint])
                   == Approx(-0.1*(1e5*(1./0.4 + 1.)*(-0.3*10 - 0.5*20) + 0.5*(10*10 + 20*20)*veloc_dot_grad_mass)));
        }
      }
    }
  }
}

TEST_CASE("derivative")
{
  const int row_size = CARTDG_MAX_BASIS_ROW_SIZE;
  cartdg::Gauss_lobatto basis (row_size);
  cartdg::Kernel_settings settings;
  Eigen::Matrix<double, row_size, row_size> diff_mat = basis.diff_mat();
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
      cartdg::variable_derivative<1, row_size>(read, write, 0, diff_mat, settings.d_pos);
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
      cartdg::variable_derivative<1, row_size>(read, write, 0, diff_mat, settings.d_pos);
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
      cartdg::variable_derivative<1, row_size>(read, write, 0, diff_mat, settings.d_pos);
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
      double coefs [] {1.103, -4.044, 0.392};
      for (int i_dim = 0; i_dim < 3; ++i_dim)
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
              for (int j_dim : {0, 1, 2}) read[i][j][k] += coefs[j_dim]*basis.node(inds[j_dim]);
            }
          }
        }
        cartdg::variable_derivative<3, row_size>(read[0][0], write[0][0], i_dim, diff_mat, settings.d_pos);
        for (int i = 0; i < row_size; ++i)
        {
          for (int j = 0; j < row_size; ++j)
          {
            for (int k = 0; k < row_size; ++k)
            {
              REQUIRE(write[i][j][k] == Approx(coefs[i_dim]));
            }
          }
        }
      }
    }
    SECTION("multi-element")
    {
      cartdg::Storage_params params {3, 3, 1, row_size};
      cartdg::elem_vec elements;
      double coefs [] {1.103, -4.044, 0.392};
      for (int i_elem = 0; i_elem < 3; ++i_elem)
      {
        elements.emplace_back(new cartdg::Element {params});
        for (int i = 0; i < row_size; ++i)
        {
          elements.back()->stage(0)[i] = coefs[i_elem]*basis.node(i);
        }
      }
      cartdg::get_local_derivative(1, row_size)(elements, 0, 0, basis, settings);
      for (int i_elem = 0; i_elem < 3; ++i_elem)
      {
        for (int i = 0; i < row_size; ++i)
        {
          REQUIRE(elements[i_elem]->derivative()[i] == Approx(coefs[i_elem]));
        }
      }
    }
    SECTION("multivariable")
    {
      cartdg::Storage_params params {3, 3, 1, row_size};
      cartdg::elem_vec elements;
      SECTION("differentiate flow vars")
      {
        for (int i_elem : {0, 1})
        {
          elements.emplace_back(new cartdg::Element {params});
          double* stage = elements[i_elem]->stage(0);
          for (int i = 0; i < row_size; ++i)
          {
            stage[i + 1*row_size] = std::pow(basis.node(i), 2);
            stage[i + 2*row_size] = std::pow(basis.node(i), 3);
          }
        }
        cartdg::get_local_derivative(1, row_size)(elements, 1, 0, basis, settings);
        for (int i_elem : {0, 1})
        {
          for (int i = 0; i < row_size; ++i)
          {
            REQUIRE(elements[i_elem]->derivative()[i] == Approx(2*basis.node(i)).margin(1e-14));
          }
        }
        cartdg::get_local_derivative(1, row_size)(elements, 2, 0, basis, settings);
        for (int i_elem : {0, 1})
        {
          for (int i = 0; i < row_size; ++i)
          {
            REQUIRE(elements[i_elem]->derivative()[i] == Approx(3*std::pow(basis.node(i), 2)).margin(1e-14));
          }
        }
      }
    }
  }
}

TEST_CASE("req_visc")
{
  const int row_size = CARTDG_MAX_BASIS_ROW_SIZE;
  const int n_qpoint = cartdg::custom_math::pow(row_size, 3);
  cartdg::Gauss_lobatto basis (row_size);
  cartdg::Kernel_settings settings;
  settings.d_pos = 0.5;
  cartdg::Storage_params params {3, 5, 3, row_size};
  cartdg::elem_vec elements;
  elements.emplace_back(new cartdg::Element {params});
  elements.emplace_back(new cartdg::Element {params});
  for (int i_elem = 0; i_elem < 2; ++i_elem)
  {
    double* stage = elements[i_elem]->stage(0);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      stage[0*n_qpoint + i_qpoint] = 0.;
      stage[1*n_qpoint + i_qpoint] = 0.;
      stage[2*n_qpoint + i_qpoint] = 0.;
      stage[3*n_qpoint + i_qpoint] = 1.225;
      stage[4*n_qpoint + i_qpoint] = 101325/0.4;
    }
  }
  elements[1]->stage(0)[3*n_qpoint + 3] = 1.5; // set an anomaly to trip the indicator in element 1
  cartdg::get_req_visc_regular_convective(3, row_size)(elements, basis, settings);
  REQUIRE(elements[0]->viscosity()[0] == 0.);
  REQUIRE(elements[0]->viscosity()[7] == 0.);
  REQUIRE(elements[1]->viscosity()[0] == Approx(0.5*340.29/(row_size - 1.)).margin(0.01));
  REQUIRE(elements[1]->viscosity()[7] == Approx(0.5*340.29/(row_size - 1.)).margin(0.01));
}

TEST_CASE("av_flux")
{
  #if CARTDG_MAX_BASIS_ROW_SIZE >= 3
  cartdg::Gauss_lobatto basis (3);
  cartdg::Kernel_settings settings;
  settings.d_pos = 0.5;
  settings.d_t_by_d_pos = 3.;
  cartdg::Storage_params params {3, 4, 2, 3};
  cartdg::elem_vec elements;
  for (int i_elem : {0, 1})
  {
    elements.emplace_back(new cartdg::Element {params});
    for (int i = 0; i < 4; ++i) elements[i_elem]->viscosity()[i] = 0.;
    for (int i = 0; i < 9; ++i) elements[i_elem]->derivative()[i] = 0.;
  }
  elements[0]->viscosity()[0] = 0.2;
  elements[1]->derivative()[1] = 0.3;
  elements[1]->derivative()[8] = 0.4;
  for (int i = 0; i < 4; ++i) elements[1]->viscosity()[i] = 1.;
  for (int i = 0; i < 9; ++i) elements[0]->derivative()[i] = 1.;
  cartdg::get_av_flux(2, 3)(elements, basis, settings);

  REQUIRE(elements[0]->derivative()[0] == Approx(0.2*1.5));
  REQUIRE(elements[0]->derivative()[1] == Approx(0.1*1.5));
  REQUIRE(elements[0]->derivative()[2] == Approx(0.0*1.5));
  REQUIRE(elements[0]->derivative()[3] == Approx(0.1*1.5));
  REQUIRE(elements[0]->derivative()[4] == Approx(0.05*1.5));
  REQUIRE(elements[0]->derivative()[5] == Approx(0.0*1.5));

  REQUIRE(elements[1]->derivative()[0] == Approx(0.0));
  REQUIRE(elements[1]->derivative()[1] == Approx(0.3*1.5));
  REQUIRE(elements[1]->derivative()[8] == Approx(0.4*1.5));
  #endif
}
