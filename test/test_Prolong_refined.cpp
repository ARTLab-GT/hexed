#include <catch2/catch_all.hpp>
#include <hexed/Spatial.hpp>
#include <hexed/Navier_stokes.hpp>
#include <hexed/Gauss_legendre.hpp>

TEST_CASE("Prolong_refined") {
  // test that prolongation operator is approximately correct for an exponential function
  const int row_size {hexed::config::max_row_size};
  hexed::Gauss_legendre basis {row_size};
  double coarse [5][row_size][row_size] {};
  double fine [2][5][row_size][row_size] {};
  std::vector<std::vector<hexed::Kernel_face_refinement>> ref_faces(1);
  hexed::Vector_view<std::vector<hexed::Kernel_face_refinement>&, std::vector<hexed::Kernel_face_refinement>> ref_face_v {ref_faces};
  for (int i_var = 0; i_var < 5; ++i_var) {
    for (int i_node = 0; i_node < row_size; ++i_node) {
      for (int j_node = 0; j_node < row_size; ++j_node) {
        coarse[i_var][i_node][j_node] = std::exp(basis.node(i_node) + 0.5*basis.node(j_node)) + i_var;
      }
    }
  }
  ref_faces.back().emplace_back();
  hexed::Kernel_face_refinement& ref = ref_faces.back().back();
  ref.coarse[0] = coarse[0][0];
  ref.coarse[1] = nullptr;
  for (int i_fine = 0; i_fine < 2; ++i_fine) {
    ref.fine[i_fine][0] = fine[i_fine][0][0];
    ref.fine[i_fine][1] = nullptr;
  }
  ref.mask = 0;

  SECTION("split dim 0") {
    ref.split_dim = 0;
    (*hexed::kernel_factory<hexed::Spatial<hexed::Navier_stokes<false>::Pde, false>::Prolong_refined>(3, row_size, basis, 0, 5))(ref_face_v);
    for (int i_half : {0, 1}) {
      for (int i_node = 0; i_node < row_size; ++i_node) {
        for (int j_node = 0; j_node < row_size; ++j_node) {
          for (int i_var = 0; i_var < 5; ++i_var) {
            double prolonged {fine[i_half][i_var][i_node][j_node]};
            double correct {std::exp((basis.node(i_node) + i_half)/2. + 0.5*basis.node(j_node)) + i_var};
            REQUIRE(prolonged == Catch::Approx(correct).margin(1e-4));
          }
        }
      }
    }
  }

  SECTION("split dim 1") {
    ref.split_dim = 1;
    (*hexed::kernel_factory<hexed::Spatial<hexed::Navier_stokes<false>::Pde, false>::Prolong_refined>(3, row_size, basis, 0, 5))(ref_face_v);
    for (int j_half : {0, 1}) {
      for (int i_node = 0; i_node < row_size; ++i_node) {
        for (int j_node = 0; j_node < row_size; ++j_node) {
          for (int i_var = 0; i_var < 5; ++i_var) {
            double prolonged {fine[j_half][i_var][i_node][j_node]};
            double correct {std::exp(basis.node(i_node) + 0.5*(basis.node(j_node) + j_half)/2.) + i_var};
            REQUIRE(prolonged == Catch::Approx(correct).margin(1e-4));
          }
        }
      }
    }
  }
}
