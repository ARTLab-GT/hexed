#include <catch2/catch.hpp>
#include <Prolong_refined.hpp>
#include <Gauss_legendre.hpp>

TEST_CASE("Prolong_refined")
{
  SECTION("no stretching")
  {
    // test that prolongation operator is approximately correct for an exponential function
    const int row_size {cartdg::config::max_row_size};
    cartdg::Storage_params params {2, 5, 3, row_size};
    cartdg::Gauss_legendre basis {row_size};
    double coarse [5][row_size][row_size] {};
    std::vector<cartdg::Refined_face> ref_faces;
    ref_faces.emplace_back(params, coarse[0][0]);
    cartdg::Vector_view<cartdg::Refined_face&, cartdg::Refined_face> ref_face_v {ref_faces};

    for (int i_var = 0; i_var < 5; ++i_var) {
      for (int i_node = 0; i_node < row_size; ++i_node) {
        for (int j_node = 0; j_node < row_size; ++j_node) {
          coarse[i_var][i_node][j_node] = std::exp(basis.node(i_node) + 0.5*basis.node(j_node)) + i_var;
        }
      }
    }

    (*cartdg::kernel_factory<cartdg::Prolong_refined>(3, row_size, basis))(ref_face_v);

    for (int i_half : {0, 1}) {
      for (int j_half : {0, 1}) {
        for (int i_node = 0; i_node < row_size; ++i_node) {
          for (int j_node = 0; j_node < row_size; ++j_node) {
            for (int i_var = 0; i_var < 5; ++i_var) {
              double prolonged {ref_faces[0].fine_face(i_half*2 + j_half)[(i_var*row_size + i_node)*row_size + j_node]};
              double correct {std::exp((basis.node(i_node) + i_half)/2. + 0.5*(basis.node(j_node) + j_half)/2.) + i_var};
              REQUIRE(prolonged == Approx(correct).margin(1e-4));
            }
          }
        }
      }
    }
  }
}
