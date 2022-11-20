#include <catch2/catch.hpp>
#include "testing_utils.hpp"
#include <hexed/Local_deformed.hpp>
#include <hexed/Gauss_legendre.hpp>

// note: this test is essentially the same as for Cartesian because anything else is pretty hard to implement.
// The real test happens at the integration level in `test_Solver`.

class Identity_basis : public hexed::Basis
{
  public:
  Identity_basis (int row_size_arg) : Basis(row_size_arg) {}
  double node(int i) const { return 0.;}
  Eigen::MatrixXd diff_mat() const
  {
    return Eigen::MatrixXd::Identity(row_size, row_size);
  }
  Eigen::VectorXd node_weights() const
  {
    return Eigen::VectorXd::Ones(row_size);
  }
  Eigen::MatrixXd boundary() const
  {
    return Eigen::MatrixXd::Zero(2, row_size);
  }
  Eigen::VectorXd orthogonal(int degree) const
  {
    return Eigen::VectorXd::Zero(row_size);
  }
};

TEST_CASE("Local_deformed")
{
  double d_t = 0.2;

  SECTION("1D")
  {
    int n_elem = 5;
    hexed::Storage_params params {2, 3, 1, 2};
    int n_qpoint = params.n_qpoint();
    std::vector<std::unique_ptr<hexed::Deformed_element>> elements;
    Identity_basis basis {int(params.row_size)};
    double mass = 1.225; double veloc = 10; double pres = 1e5;
    double mmtm = mass*veloc; double ener = pres/0.4 + 0.5*mass*veloc*veloc;
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new hexed::Deformed_element {params, {}, 2.});
      double* state = elements[i_elem]->stage(0);
      elements[i_elem]->set_jacobian(hexed::Gauss_legendre(params.row_size));
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        state[0*n_qpoint + i_qpoint] = mmtm;
        state[1*n_qpoint + i_qpoint] = mass;
        state[2*n_qpoint + i_qpoint] = ener;
        // set RK reference state
        for (int i_var = 0; i_var < 3; ++i_var) state[(3 + i_var)*n_qpoint + i_qpoint] = 30.;
      }
      for (int i_face = 0; i_face < n_qpoint/params.row_size*3*2*1; ++i_face) {
        elements[i_elem]->face()[i_face] = 0.;
      }
    }
    def_elem_view elem_view {elements};
    (*hexed::kernel_factory<hexed::Local_deformed>(1, 2, basis, d_t*.9, .9, .1))(elem_view);
    for (auto& element : elements)
    {
      double* state = element->stage(0);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        REQUIRE(state[0*n_qpoint + i_qpoint] == Approx(0.1*30. + 0.9*(mmtm - 0.1*(mmtm*mmtm/mass + pres))));
        REQUIRE(state[1*n_qpoint + i_qpoint] == Approx(0.1*30. + 0.9*(mass - 0.1*mmtm)));
        REQUIRE(state[2*n_qpoint + i_qpoint] == Approx(0.1*30. + 0.9*(ener - 0.1*((ener + pres)*veloc))));
      }
    }
  }

  SECTION("3D")
  {
    const int n_elem = 5;
    hexed::Storage_params params {2, 5, 3, 3};
    int n_qpoint = params.n_qpoint();
    std::vector<std::unique_ptr<hexed::Deformed_element>> elements;
    Identity_basis basis (params.row_size);
    double mass = 1.225;
    double veloc0 = 10; double veloc1 = -20; double veloc2 = 30;
    double pres = 1e5;
    double ener = pres/0.4 + 0.5*mass*(veloc0*veloc0 + veloc1*veloc1 + veloc2*veloc2);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new hexed::Deformed_element {params, {}, 2.});
      double* state = elements[i_elem]->stage(0);
      elements[i_elem]->set_jacobian(hexed::Gauss_legendre(params.row_size)); // just sets jacobian to identity
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        state[0*n_qpoint + i_qpoint] = mass*veloc0;
        state[1*n_qpoint + i_qpoint] = mass*veloc1;
        state[2*n_qpoint + i_qpoint] = mass*veloc2;
        state[3*n_qpoint + i_qpoint] = mass;
        state[4*n_qpoint + i_qpoint] = ener;
        // set RK reference state (since RK weight is 1., inializing to 0 works fine)
        for (int i_var = 0; i_var < 5; ++i_var) state[(5 + i_var)*n_qpoint + i_qpoint] = 0.;
        // set Jacobian to identity
      }
      for (int i_face = 0; i_face < n_qpoint/params.row_size*5*2*3; ++i_face) {
        elements[i_elem]->face()[i_face] = 0.;
      }
      elements[i_elem]->time_step_scale()[1] = 0.16;
    }
    def_elem_view elem_view {elements};
    (*hexed::kernel_factory<hexed::Local_deformed>(3, params.row_size, basis, d_t, 1., 0.))(elem_view);
    for (auto& element : elements)
    {
      double* state = element->stage(0);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        double tss = (i_qpoint == 1) ? 0.16 : 1.;
        REQUIRE(state[0*n_qpoint + i_qpoint] - mass*veloc0 == Approx(-0.1*tss*(mass*veloc0*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE(state[1*n_qpoint + i_qpoint] - mass*veloc1 == Approx(-0.1*tss*(mass*veloc1*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE(state[2*n_qpoint + i_qpoint] - mass*veloc2 == Approx(-0.1*tss*(mass*veloc2*(veloc0 + veloc1 + veloc2) + pres)));
        REQUIRE(state[3*n_qpoint + i_qpoint] - mass        == Approx(-0.1*tss*mass*(veloc0 + veloc1 + veloc2)));
        REQUIRE(state[4*n_qpoint + i_qpoint] - ener        == Approx(-0.1*tss*((ener + pres)*(veloc0 + veloc1 + veloc2))));
      }
    }
  }

  SECTION("2D non-constant")
  {
    const int n_elem = 5;
    const int row_size = std::min<int>(6, int(hexed::config::max_row_size));
    hexed::Storage_params params {3, 4, 2, row_size};
    int n_qpoint = params.n_qpoint();
    std::vector<std::unique_ptr<hexed::Deformed_element>> elements;
    hexed::Gauss_legendre basis (row_size);

    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      elements.emplace_back(new hexed::Deformed_element {params, {}, 2.});
      double* state = elements[i_elem]->stage(0);
      double* face = elements[i_elem]->face();
      elements[i_elem]->set_jacobian(hexed::Gauss_legendre(params.row_size));

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
          // set initial state
          state[0*n_qpoint + i_qpoint] = mass*veloc[0];
          state[1*n_qpoint + i_qpoint] = mass*veloc[1];
          state[2*n_qpoint + i_qpoint] = mass;
          state[3*n_qpoint + i_qpoint] = ener;
          // set RK reference state to negative of initial state
          for (int i_var = 0; i_var < 4; ++i_var) {
            state[(4 + i_var)*n_qpoint + i_qpoint] = -state[i_var*n_qpoint + i_qpoint];
          }
        }
        // set face state to match interior state
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
    def_elem_view elem_view {elements};
    (*hexed::kernel_factory<hexed::Local_deformed>(2, row_size, basis, d_t, 0, 0))(elem_view);
    for (auto& element : elements) {
      double* state = element->stage(0);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        REQUIRE(state[2*n_qpoint + i_qpoint] == Approx(-0.1*(0.1*10 - 0.2*20)));
        REQUIRE(state[3*n_qpoint + i_qpoint] == Approx(-0.1*(1e5*(1./0.4 + 1.)*(-0.3*10 - 0.5*20)
                                                       + 0.5*(10*10 + 20*20)*(0.1*10 - 0.2*20))));
      }
    }
  }
}
