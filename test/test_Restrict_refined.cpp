#include <catch2/catch.hpp>
#include <hexed/Restrict_refined.hpp>
#include <hexed/Gauss_legendre.hpp>

TEST_CASE("Restrict_refined")
{
  // test that restriction operator is approximately correct for an exponential function
  const int row_size {hexed::config::max_row_size};
  hexed::Storage_params params {2, 5, 3, row_size};
  hexed::Gauss_legendre basis {row_size};
  double coarse [5][row_size][row_size] {};
  std::vector<hexed::Refined_face> ref_faces;
  hexed::Vector_view<hexed::Refined_face&, hexed::Refined_face> ref_face_v {ref_faces};
  auto check = [&](double factor) {
    for (int i_var = 0; i_var < 5; ++i_var) {
      for (int i_node = 0; i_node < row_size; ++i_node) {
        for (int j_node = 0; j_node < row_size; ++j_node) {
          double restricted {coarse[i_var][i_node][j_node]};
          double correct {std::exp(basis.node(i_node) + 0.5*basis.node(j_node)) + i_var};
          REQUIRE(restricted == Approx(factor*correct).margin(1e-4));
        }
      }
    }
  };

  SECTION("no stretching")
  {
    ref_faces.emplace_back(params, coarse[0][0]);
    for (int i_half : {0, 1}) {
      for (int j_half : {0, 1}) {
        for (int i_node = 0; i_node < row_size; ++i_node) {
          for (int j_node = 0; j_node < row_size; ++j_node) {
            for (int i_var = 0; i_var < 5; ++i_var) {
              ref_faces[0].fine_face(i_half*2 + j_half)[(i_var*row_size + i_node)*row_size + j_node]
                = std::exp((basis.node(i_node) + i_half)/2. + 0.5*(basis.node(j_node) + j_half)/2.) + i_var;
            }
          }
        }
      }
    }
    (*hexed::kernel_factory<hexed::Restrict_refined>(3, row_size, basis))(ref_face_v);
    check(1.);
  }

  SECTION("stretch dimension 0")
  {
    ref_faces.emplace_back(params, coarse[0][0], std::array<bool, 2>{true, false});
    for (int j_half : {0, 1}) {
      for (int i_node = 0; i_node < row_size; ++i_node) {
        for (int j_node = 0; j_node < row_size; ++j_node) {
          for (int i_var = 0; i_var < 5; ++i_var) {
            ref_faces[0].fine_face(j_half)[(i_var*row_size + i_node)*row_size + j_node]
              = std::exp(basis.node(i_node) + 0.5*(basis.node(j_node) + j_half)/2.) + i_var;
          }
        }
      }
    }
    (*hexed::kernel_factory<hexed::Restrict_refined>(3, row_size, basis))(ref_face_v);
    check(.5);
  }

  SECTION("stretch dimension 1")
  {
    ref_faces.emplace_back(params, coarse[0][0], std::array<bool, 2>{false, true});
    for (int i_half : {0, 1}) {
      for (int i_node = 0; i_node < row_size; ++i_node) {
        for (int j_node = 0; j_node < row_size; ++j_node) {
          for (int i_var = 0; i_var < 5; ++i_var) {
            ref_faces[0].fine_face(i_half)[(i_var*row_size + i_node)*row_size + j_node]
              = std::exp((basis.node(i_node) + i_half)/2. + 0.5*basis.node(j_node)) + i_var;
          }
        }
      }
    }
    (*hexed::kernel_factory<hexed::Restrict_refined>(3, row_size, basis))(ref_face_v);
    check(.5);
  }
}
