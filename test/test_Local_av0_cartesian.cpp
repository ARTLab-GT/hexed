#include <catch2/catch.hpp>
#include "testing_utils.hpp"
#include <hexed/Local_av0_cartesian.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/Qpoint_func.hpp>
#include <hexed/Write_face.hpp>

TEST_CASE("Local_av0_cartesian.hpp")
{
  const int row_size = hexed::config::max_row_size;
  static_assert(row_size >= 5); // must be able to represent quartic polynomial exactly
  hexed::Storage_params params {2, 5, 3, row_size};
  const int n_qpoint = params.n_qpoint();
  hexed::Gauss_legendre basis(row_size);
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
  std::vector<std::unique_ptr<hexed::Element>> elems;
  elems.emplace_back(new hexed::Element(params, {}, .4, 1));
  auto& elem = *elems.back();
  double scale [] {-.2, .3, .7};
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
  car_elem_view elem_view {elems};
  (hexed::Write_face<3, row_size>(basis))(elem_view); // `Local_av0_cartesian` expects faces to contain numerical LDG state
  (hexed::Local_av0_cartesian<3, row_size>(basis, .3, 1.))(elem_view);
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
