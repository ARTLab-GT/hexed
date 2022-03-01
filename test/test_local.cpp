#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <get_write_face.hpp>
#include <get_local_convective.hpp>
#include <get_req_visc_regular_convective.hpp>
#include <get_local_gradient.hpp>
#include <Gauss_lobatto.hpp>
#include <Gauss_legendre.hpp>
#include <Equidistant.hpp>
#include <math.hpp>
#include <Kernel_settings.hpp>

class Identity_basis : public cartdg::Basis
{
  public:
  Identity_basis (int row_size_arg) : Basis(row_size_arg) {}
  double node(int i) { return 0.;}
  Eigen::MatrixXd diff_mat()
  {
    return Eigen::MatrixXd::Identity(row_size, row_size);
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

TEST_CASE("write_face")
{
  if (cartdg::config::max_row_size >= 5)
  {
    const int row_size {5};
    cartdg::Kernel_settings settings;
    cartdg::Equidistant basis {row_size};
    cartdg::Storage_params params {2, 5, 3, row_size};
    cartdg::elem_vec elements;
    elements.emplace_back(new cartdg::Element {params});
    const int n_qpoint {params.n_qpoint()};
    for (int i_var : {0, 1})
    {
      for (int i_row = 0; i_row < row_size; ++i_row)
      for (int j_row = 0; j_row < row_size; ++j_row)
      for (int k_row = 0; k_row < row_size; ++k_row)
      {
        int i_qpoint = (i_row*row_size + j_row)*row_size + k_row;
        double value = 0.1*i_var + 0.2*basis.node(i_row) + 0.3*basis.node(j_row) + 0.4*basis.node(k_row);
        elements[0]->stage(settings.i_read)[i_var*n_qpoint + i_qpoint] = value;
      }
    }
    cartdg::get_write_face(3, row_size)(elements, basis, settings);
    double* face = elements[0]->face();
    const int n_face {params.n_dof()/row_size};
    REQUIRE(face[0*n_face + 0*n_qpoint/row_size + 0] == Approx(0.).margin(1e-10));
    REQUIRE(face[0*n_face + 0*n_qpoint/row_size + 3] == Approx(0.75*0.4));
    REQUIRE(face[0*n_face + 0*n_qpoint/row_size + 5] == Approx(0.25*0.3));
    REQUIRE(face[0*n_face + 0*n_qpoint/row_size + 6] == Approx(0.25*(0.3 + 0.4)));
    REQUIRE(face[0*n_face + 1*n_qpoint/row_size + 0] == Approx(0.1));
    REQUIRE(face[1*n_face + 0*n_qpoint/row_size + 0] == Approx(1.*0.2));
    REQUIRE(face[1*n_face + 0*n_qpoint/row_size + 1] == Approx(1.*0.2 + 0.25*0.4));
    REQUIRE(face[2*n_face + 0*n_qpoint/row_size + 0] == Approx(0.).margin(1e-10));
    REQUIRE(face[2*n_face + 0*n_qpoint/row_size + 1] == Approx(0.25*0.4));
    REQUIRE(face[2*n_face + 0*n_qpoint/row_size + 5] == Approx(0.25*0.2));
    REQUIRE(face[5*n_face + 0*n_qpoint/row_size + 8] == Approx(1.*0.4 + 0.25*0.2 + 0.75*0.3));
  }
}

TEST_CASE("Local convective")
{
  cartdg::Kernel_settings settings;
  settings.d_t_by_d_pos = 0.1;
  SECTION("1D")
  {
    int n_elem = 5;
    cartdg::Storage_params params {2, 3, 1, 2};
    int n_qpoint = params.n_qpoint();
    cartdg::elem_vec elements;
    Identity_basis basis {int(params.row_size)};
    double mass = 1.225; double veloc = 10; double pres = 1e5;
    double mmtm = mass*veloc; double ener = pres/0.4 + 0.5*mass*veloc*veloc;
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new cartdg::Element {params});
      double* read = elements[i_elem]->stage(settings.i_read);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
          read[0*n_qpoint + i_qpoint] = mmtm;
          read[1*n_qpoint + i_qpoint] = mass;
          read[2*n_qpoint + i_qpoint] = ener;
      }
      for (int i_face = 0; i_face < n_qpoint/params.row_size*3*2*1; ++i_face)
      {
        elements[i_elem]->face()[i_face] = 0.;
      }
    }
    cartdg::get_local_convective(1, 2)(elements, basis, settings);
    for (auto& element : elements)
    {
      double* r = element->stage(settings.i_read);
      double* w = element->stage(settings.i_write);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        REQUIRE(w[0*n_qpoint + i_qpoint] - r[0*n_qpoint + i_qpoint] == Approx(-0.1*(mmtm*mmtm/mass + pres)));
        REQUIRE(w[1*n_qpoint + i_qpoint] - r[1*n_qpoint + i_qpoint] == Approx(-0.1*mmtm));
        REQUIRE(w[2*n_qpoint + i_qpoint] - r[2*n_qpoint + i_qpoint] == Approx(-0.1*((ener + pres)*veloc)));
      }
    }
  }

  SECTION("3D")
  {
    const int n_elem = 5;
    cartdg::Storage_params params {2, 5, 3, 3};
    int n_qpoint = params.n_qpoint();
    cartdg::elem_vec elements;
    Identity_basis basis (params.row_size);
    double mass = 1.225;
    double veloc0 = 10; double veloc1 = -20; double veloc2 = 30;
    double pres = 1e5;
    double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1 + veloc2*veloc2);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new cartdg::Element {params});
      double* read = elements[i_elem]->stage(settings.i_read);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
          read[0*n_qpoint + i_qpoint] = mass*veloc0;
          read[1*n_qpoint + i_qpoint] = mass*veloc1;
          read[2*n_qpoint + i_qpoint] = mass*veloc2;
          read[3*n_qpoint + i_qpoint] = mass;
          read[4*n_qpoint + i_qpoint] = ener;
      }
      for (int i_face = 0; i_face < n_qpoint/params.row_size*5*2*3; ++i_face)
      {
        elements[i_elem]->face()[i_face] = 0.;
      }
      elements[i_elem]->time_step_scale()[1] = 0.16;
    }
    cartdg::get_local_convective(3, 3)(elements, basis, settings);
    for (auto& element : elements)
    {
      double* r = element->stage(settings.i_read);
      double* w = element->stage(settings.i_write);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        double tss = (i_qpoint == 1) ? 0.16 : 1.;
        REQUIRE(w[0*n_qpoint + i_qpoint] - r[0*n_qpoint + i_qpoint]
                == Approx(-0.1*tss*(mass*veloc0*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE(w[1*n_qpoint + i_qpoint] - r[1*n_qpoint + i_qpoint]
                == Approx(-0.1*tss*(mass*veloc1*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE(w[2*n_qpoint + i_qpoint] - r[2*n_qpoint + i_qpoint]
                == Approx(-0.1*tss*(mass*veloc2*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE(w[3*n_qpoint + i_qpoint] - r[3*n_qpoint + i_qpoint]
                == Approx(-0.1*tss*mass*(veloc0 + veloc1 + veloc2)));
        REQUIRE(w[4*n_qpoint + i_qpoint] - r[4*n_qpoint + i_qpoint]
                == Approx(-0.1*tss*((ener + pres)*(veloc0 + veloc1 + veloc2))));
      }
    }
  }

  SECTION("2D non-constant")
  {
    const int n_elem = 5;
    const int row_size = std::min<int>(6, int(cartdg::config::max_row_size));
    cartdg::Storage_params params {3, 4, 2, row_size};
    int n_qpoint = params.n_qpoint();
    cartdg::elem_vec elements;
    cartdg::Gauss_legendre basis (row_size);

    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new cartdg::Element {params});
      double* read = elements[i_elem]->stage(settings.i_read);
      double* face = elements[i_elem]->face();

      #define SET_VARS \
        double mass = 1 + 0.1*pos0 + 0.2*pos1; \
        double veloc [] {10, -20}; \
        double pres = 1e5*(1. - 0.3*pos0 + 0.5*pos1); \
        double ener = pres/0.4 + 0.5*mass*(veloc[0]*veloc[0] + veloc[1]*veloc[1]); \

      for (int i = 0; i < row_size; ++i)
      {
        for (int j = 0; j < row_size; ++j)
        {
          int i_qpoint = i*row_size + j;
          double pos0 = basis.node(i); double pos1 = basis.node(j);
          SET_VARS
          read[0*n_qpoint + i_qpoint] = mass*veloc[0];
          read[1*n_qpoint + i_qpoint] = mass*veloc[1];
          read[2*n_qpoint + i_qpoint] = mass;
          read[3*n_qpoint + i_qpoint] = ener;
        }
        for (int i_dim : {0, 1})
        {
          for (int positive : {0, 1})
          {
            double* qpoint_start = face + (2*i_dim + positive)*row_size*4 + i;
            double pos0 {i_dim ? basis.node(i) : positive};
            double pos1 {i_dim ? positive : basis.node(i)};
            SET_VARS
            qpoint_start[0*row_size] = mass*veloc[0]*veloc[i_dim];
            qpoint_start[1*row_size] = mass*veloc[1]*veloc[i_dim];
            qpoint_start[2*row_size] = mass*veloc[i_dim];
            qpoint_start[3*row_size] = (ener + pres)*veloc[i_dim];
            qpoint_start[i_dim*row_size] += pres;
          }
        }
      }
      #undef SET_VARS
    }
    cartdg::get_local_convective(2, row_size)(elements, basis, settings);
    for (auto& element : elements)
    {
      double* r = element->stage(settings.i_read);
      double* w = element->stage(settings.i_write);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
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

TEST_CASE("req_visc")
{
  const int row_size = cartdg::config::max_row_size;
  const int n_qpoint = cartdg::custom_math::pow(row_size, 3);
  cartdg::Gauss_legendre basis (row_size);
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

TEST_CASE("local_gradient")
{
  const int row_size = cartdg::config::max_row_size;
  static_assert (row_size >= 3); // this test requires at least quadratic polynomial
  cartdg::Storage_params params {3, 5, 3, row_size};
  cartdg::elem_vec elements;
  elements.emplace_back(new cartdg::Element {params});
  cartdg::Element& element = *elements.back();
  const int n_qpoint = cartdg::custom_math::pow(row_size, 3);
  cartdg::Gauss_legendre basis (row_size);
  double nodes [n_qpoint][3];
  // compute the positions of the 3D quadrature points
  for (int i_row = 0; i_row < row_size; ++i_row) {
    for (int j_row = 0; j_row < row_size; ++j_row) {
      for (int k_row = 0; k_row < row_size; ++k_row) {
        int i_qpoint = (i_row*row_size + j_row)*row_size + k_row;
        int row [3] {i_row, j_row, k_row};
        for (int i_dim = 0; i_dim < 3; ++i_dim) nodes[i_qpoint][i_dim] = basis.node(row[i_dim]);
      }
    }
  }
  // write an arbitrary triquadratic polynomial to the mass
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    double* node = nodes[i_qpoint];
    element.stage(0)[3*n_qpoint + i_qpoint] = 2.*node[0]*node[0] - 3.*node[0]*node[1]*node[1] + 1.*node[0]*node[1]*node[2];
  }
  cartdg::Kernel_settings settings;
  cartdg::get_local_gradient(3, row_size)(elements, 3, basis, settings);
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    double* node = nodes[i_qpoint];
    // check that derivatives match analytic
    REQUIRE(element.stage(1)[0*n_qpoint + i_qpoint] == Approx(4.*node[0] - 3.*node[1]*node[1] + 1.*node[1]*node[2]));
    REQUIRE(element.stage(1)[1*n_qpoint + i_qpoint] == Approx(-6.*node[0]*node[1] + 1.*node[0]*node[2]));
    REQUIRE(element.stage(1)[2*n_qpoint + i_qpoint] == Approx(1.*node[0]*node[1]));
  }
}
