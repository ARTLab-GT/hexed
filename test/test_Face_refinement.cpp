#include <catch2/catch_all.hpp>
#include <hexed/Face_refinement.hpp>

TEST_CASE("Face_refinement") {
  hexed::Storage_params params {2, 5, 3, 2};
  auto f = std::make_unique<hexed::Face>(params, 0, 1);
  REQUIRE_THROWS(hexed::Face_refinement(*f, 2)); // split_dim is too large
  hexed::Face_refinement fr(*f, 1);
  REQUIRE(&fr.coarse() == f.get());
  REQUIRE(fr.split_dim() == 1);
  REQUIRE(fr.alive());
  REQUIRE(fr.fine()[0]->i_dim() == f->i_dim());
  REQUIRE(fr.fine()[1]->sign() == f->sign());
  REQUIRE(f->connected());
  REQUIRE(!f->associated());
  REQUIRE(f->face_ref_fine() == &fr);
  REQUIRE(f->face_ref_coarse() == nullptr);
  REQUIRE(fr.fine()[0]->associated());
  REQUIRE(!fr.fine()[0]->connected());
  REQUIRE(fr.fine()[1]->face_ref_fine() == nullptr);
  REQUIRE(fr.fine()[1]->face_ref_coarse() == &fr);
  f.reset();
  REQUIRE(!fr.alive());
  REQUIRE_THROWS(fr.coarse());
}
