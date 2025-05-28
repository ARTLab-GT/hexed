#include <catch2/catch_all.hpp>
#include <hexed/Face_refinement.hpp>
#include <hexed/vertex_inds.hpp>

TEST_CASE("Face_refinement") {
  hexed::Storage_params params {2, 5, 3, 2};
  auto f = std::make_unique<hexed::Face>(params, 0, 1, false);
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
  REQUIRE(fr.is_deformed() == false);
  f.reset();
  REQUIRE(!fr.alive());
  REQUIRE_THROWS(fr.coarse());
  REQUIRE(fr.is_deformed() == false);

  SECTION("elements") {
    hexed::Storage_params params {2, 5, 3, 2};
    std::vector<std::unique_ptr<hexed::Element>> elems;
    std::vector<hexed::Face_refinement> face_refs;
    std::vector<hexed::Neighbor_connection> neighb_cons;
    SECTION("2 on 2") {
      for (int i = 0; i < 4; ++i) elems.push_back(std::make_unique<hexed::Element>(params));
      face_refs.emplace_back(elems[0]->face(1), 0);
      face_refs.emplace_back(elems[1]->face(1), 0);
      face_refs.emplace_back(elems[2]->face(0), 1);
      face_refs.emplace_back(elems[3]->face(0), 1);
      hexed::Connection_direction dir {{0, 0}, {1, 0}, 2};
      auto inds = hexed::face_vertex_inds(3, dir);
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          int ind = inds[2*i + j];
          std::array<hexed::Face*, 2> face_arr {face_refs[i].fine()[j], face_refs[2 + ind%2].fine()[ind/2]};
          neighb_cons.emplace_back(params, face_arr, 2);
        }
      }
      for (int i_ref = 0; i_ref < (int)face_refs.size(); ++i_ref) {
        auto ref_elems = face_refs[i_ref].elements();
        REQUIRE_THAT(ref_elems[0], Catch::Matchers::RangeEquals(std::vector<hexed::Element*> {
          elems[0].get(), elems[1].get(), elems[0].get(), elems[1].get(),
        }));
        REQUIRE_THAT(ref_elems[1], Catch::Matchers::RangeEquals(std::vector<hexed::Element*> {
          elems[2].get(), elems[2].get(), elems[3].get(), elems[3].get(),
        }));
        REQUIRE(face_refs[i_ref].get_direction() == hexed::Connection_direction{{0, 0}, {1, 0}, 2});
      }
    }
    SECTION("1 on 4") {
      for (int i = 0; i < 5; ++i) elems.push_back(std::make_unique<hexed::Element>(params));
      face_refs.emplace_back(elems[4]->face(4), 0);
      face_refs.emplace_back(*face_refs[0].fine()[0], 1);
      face_refs.emplace_back(*face_refs[0].fine()[1], 1);
      hexed::Connection_direction dir {{2, 1}, {0, 1}};
      auto inds = hexed::face_vertex_inds(3, dir);
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          std::array<hexed::Face*, 2> face_arr {face_refs[1 + i].fine()[j], &elems[inds[2*i + j]]->face(3)};
          neighb_cons.emplace_back(params, face_arr);
        }
      }
      for (int i_ref = 0; i_ref < (int)face_refs.size(); ++i_ref) {
        auto ref_elems = face_refs[i_ref].elements();
        REQUIRE_THAT(ref_elems[0], Catch::Matchers::RangeEquals(std::vector<hexed::Element*>(4, elems[4].get())));
        REQUIRE_THAT(ref_elems[1], Catch::Matchers::RangeEquals(std::vector<hexed::Element*> {
          elems[0].get(), elems[1].get(), elems[2].get(), elems[3].get(),
        }));
        REQUIRE(face_refs[i_ref].get_direction() == dir);
      }
    }
    SECTION("4 on 1") {
      for (int i = 0; i < 5; ++i) elems.push_back(std::make_unique<hexed::Element>(params));
      face_refs.emplace_back(elems[4]->face(5), 0);
      face_refs.emplace_back(*face_refs[0].fine()[0], 1);
      face_refs.emplace_back(*face_refs[0].fine()[1], 1);
      hexed::Connection_direction dir {{0, 2}, {1, 1}};
      auto inds = hexed::face_vertex_inds(3, dir);
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          int ind = inds[2*i + j];
          std::array<hexed::Face*, 2> face_arr {&elems[2*i + j]->face(1), face_refs[1 + ind/2].fine()[ind%2]};
          neighb_cons.emplace_back(params, face_arr);
        }
      }
      for (int i_ref = 0; i_ref < (int)face_refs.size(); ++i_ref) {
        auto ref_elems = face_refs[i_ref].elements();
        REQUIRE_THAT(ref_elems[0], Catch::Matchers::RangeEquals(std::vector<hexed::Element*> {
          elems[0].get(), elems[1].get(), elems[2].get(), elems[3].get(),
        }));
        REQUIRE_THAT(ref_elems[1], Catch::Matchers::RangeEquals(std::vector<hexed::Element*>(4, elems[4].get())));
        REQUIRE(face_refs[i_ref].get_direction() == dir);
      }
    }
  }
}
