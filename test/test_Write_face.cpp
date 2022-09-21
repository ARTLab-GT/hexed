#include <catch2/catch.hpp>
#include <hexed/config.hpp>
#include <hexed/Write_face.hpp>
#include <hexed/Equidistant.hpp>
#include "testing_utils.hpp"

TEST_CASE("Write_face")
{
  static_assert (hexed::config::max_row_size >= 5); // this test was written for row size 5
  const int row_size {5};
  hexed::Equidistant basis {row_size};
  hexed::Storage_params params {2, 5, 3, row_size};
  std::vector<std::unique_ptr<hexed::Element>> elements;
  elements.emplace_back(new hexed::Element {params});
  const int n_qpoint {params.n_qpoint()};
  for (int i_var : {0, 1})
  {
    for (int i_row = 0; i_row < row_size; ++i_row)
    for (int j_row = 0; j_row < row_size; ++j_row)
    for (int k_row = 0; k_row < row_size; ++k_row)
    {
      int i_qpoint = (i_row*row_size + j_row)*row_size + k_row;
      double value = 0.1*i_var + 0.2*basis.node(i_row) + 0.3*basis.node(j_row) + 0.4*basis.node(k_row);
      elements[0]->stage(0)[i_var*n_qpoint + i_qpoint] = value;
    }
  }
  car_elem_view elem_view {elements};
  (*hexed::kernel_factory<hexed::Write_face>(3, row_size, basis))(elem_view);
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
