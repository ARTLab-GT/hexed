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

  SECTION("elements") {
    hexed::Storage_params params {2, 5, 3, 2};
    std::vector<std::unique_ptr<hexed::Element>> elems;
    for (int i = 0; i < 4; ++i) elems.push_back(std::make_unique<hexed::Element>(params));
    std::vector<hexed::Face_refinement> face_refs;
    face_refs.emplace_back(elems[0]->face(1), 0);
    face_refs.emplace_back(elems[1]->face(1), 0);
    face_refs.emplace_back(elems[2]->face(0), 1);
    face_refs.emplace_back(elems[3]->face(0), 1);
    std::vector<hexed::Neighbor_connection> neighb_cons;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        std::array<hexed::Face*, 2> face_arr {face_refs[i].fine()[j], face_refs[2 + j].fine()[i]};
        neighb_cons.emplace_back(params, face_arr);
      }
    }
    for (int i_ref = 0; i_ref < (int)face_refs.size(); ++i_ref) {
      auto ref_elems = face_refs[i_ref].elements();
      REQUIRE_THAT(ref_elems[i_ref >= 2], Catch::Matchers::RangeEquals(std::vector<hexed::Element*> {
        elems[0].get(), elems[1].get(), elems[0].get(), elems[1].get(),
      }));
      REQUIRE_THAT(ref_elems[i_ref < 2], Catch::Matchers::RangeEquals(std::vector<hexed::Element*> {
        elems[2].get(), elems[2].get(), elems[3].get(), elems[3].get(),
      }));
    }
  }
}
