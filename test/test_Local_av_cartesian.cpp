#include <catch2/catch.hpp>
#include "testing_utils.hpp"
#include <hexed/Local_av0_cartesian.hpp>
#include <hexed/Local_av1_cartesian.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/Qpoint_func.hpp>
#include <hexed/Write_face.hpp>

TEST_CASE("local cartesian artificial viscosity")
{
  const int row_size = hexed::config::max_row_size;
  static_assert(row_size >= 5); // must be able to represent quartic polynomial exactly
  hexed::Storage_params params {2, 5, 3, row_size};
  const int n_qpoint = params.n_qpoint();
  hexed::Gauss_legendre basis(row_size);
  std::vector<std::unique_ptr<hexed::Element>> elems;
  elems.emplace_back(new hexed::Element(params, {}, .4, 1));
  auto& elem = *elems.back();
  Eigen::VectorXd linear(row_size);
  Eigen::VectorXd cubic(row_size);
  Eigen::VectorXd deriv_cubic(row_size);
  Eigen::VectorXd d2_cubic(row_size);
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    double pos = basis.node(i_qpoint);
    linear(i_qpoint) = pos;
    cubic(i_qpoint) = pos*pos/2. - pos*pos*pos/3.;
    deriv_cubic(i_qpoint) = pos - pos*pos;
    d2_cubic(i_qpoint) = 1. - 2.*pos;
  }
  double scale [] {-.2, .3, .7};
  car_elem_view elem_view {elems};

  SECTION("interior term")
  {
    for (int i = 0; i < row_size; ++i) {
      for (int j = 0; j < row_size; ++j) {
        for (int k = 0; k < row_size; ++k) {
          int inds [] {i, j, k};
          double value = 0.;
          for (int i_dim = 0; i_dim < 3; ++i_dim) value += scale[i_dim]*cubic(inds[i_dim]);
          int i_qpoint = (i*row_size + j)*row_size + k;
          for (int i_var = 0; i_var < 5; ++i_var) {
            for (int i_stage = 0; i_stage < 2; ++i_stage) {
              elem.stage(i_stage)[i_var*n_qpoint + i_qpoint] = 0.;
            }
          }
          for (int i_stage = 0; i_stage < 2; ++i_stage) {
            elem.stage(i_stage)[1*n_qpoint + i_qpoint] = value;
          }
          elem.art_visc_coef()[i_qpoint] = linear(k);
        }
      }
    }
    (hexed::Write_face<3, row_size>(basis))(elem_view); // `Local_av0_cartesian` expects faces to contain numerical LDG state
    (hexed::Local_av0_cartesian<3, row_size>(basis, .3))(elem_view);
    for (int i = 0; i < row_size; ++i) {
      for (int j = 0; j < row_size; ++j) {
        for (int k = 0; k < row_size; ++k) {
          int inds [] {i, j, k};
          double correct = 0.;
          for (int i_dim = 0; i_dim < 3; ++i_dim) correct += linear(k)*scale[i_dim]*d2_cubic(inds[i_dim]);
          correct += scale[2]*deriv_cubic(k);
          correct *= .3/(.4*.4/2/2); // .3 for time step, .4 for mesh size, 2 for ref level
          int i_qpoint = (i*row_size + j)*row_size + k;
          REQUIRE(hexed::Physical_update()(elem, basis, i_qpoint, 0)[1] == Approx(correct).scale(1.));
        }
      }
    }
  }

  SECTION("face terms")
  {
    for (int i = 0; i < row_size; ++i) {
      for (int j = 0; j < row_size; ++j) {
        for (int k = 0; k < row_size; ++k) {
          int inds [] {i, j, k};
          double value = 0.;
          for (int i_dim = 0; i_dim < 3; ++i_dim) value += scale[i_dim]*linear(inds[i_dim]);
          int i_qpoint = (i*row_size + j)*row_size + k;
          for (int i_var = 0; i_var < 5; ++i_var) {
            for (int i_stage = 0; i_stage < 2; ++i_stage) {
              elem.stage(i_stage)[i_var*n_qpoint + i_qpoint] = 0.;
            }
          }
          for (int i_stage = 0; i_stage < 2; ++i_stage) {
            elem.stage(i_stage)[i_qpoint] = value;
          }
          elem.art_visc_coef()[i_qpoint] = .73;
        }
      }
    }
    (hexed::Write_face<3, row_size>(basis))(elem_view); // `Local_av0_cartesian` expects faces to contain numerical LDG state
    double face_data [6*5*n_qpoint/row_size];
    for (int i_data = 0; i_data < 6*5*n_qpoint/row_size; ++i_data) {
      face_data[i_data] = elem.face()[i_data];
    }
    (hexed::Local_av0_cartesian<3, row_size>(basis, .314))(elem_view);
    SECTION("flux writing")
    {
      int i_row = row_size/2; // set an arbitrary quadrature point to sample on each face
      REQUIRE(elem.face()[0*5*n_qpoint/row_size + i_row*(row_size + 1)] == Approx(.73/.4*2*-scale[0]).scale(1.));
      REQUIRE(elem.face()[1*5*n_qpoint/row_size + i_row*(row_size + 1)] == Approx(.73/.4*2*-scale[0]).scale(1.));
      REQUIRE(elem.face()[2*5*n_qpoint/row_size + i_row*(row_size + 1)] == Approx(.73/.4*2*-scale[1]).scale(1.));
      REQUIRE(elem.face()[4*5*n_qpoint/row_size + i_row*(row_size + 1)] == Approx(.73/.4*2*-scale[2]).scale(1.));
    }
    SECTION("face local term")
    {
      (hexed::Local_av1_cartesian<3, row_size>(basis, .314))(elem_view);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        // correct answer is 0 because 2nd derivative is 0
        // if face term is not correctly added answer may be nonzero
        REQUIRE(hexed::Physical_update()(elem, basis, i_qpoint, 0)[0] == Approx(0.).scale(1.));
      }
      // check that the correct state has been written to the faces
      for (int i_data = 0; i_data < 6*5*n_qpoint/row_size; ++i_data) {
        REQUIRE(elem.face()[i_data] == Approx(face_data[i_data]).scale(1.));
      }
    }
  }
}