#include <catch2/catch_all.hpp>
#include <hexed/config.hpp>
#include <hexed/Accessible_mesh.hpp>
#include <hexed/Simplex_geom.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/Printer.hpp>

TEST_CASE("Tree meshing", "[.slow]") {
  hexed::Accessible_mesh mesh({1, 5, 3, hexed::config::max_row_size}, .7, hexed::laminar);
  REQUIRE_THROWS(mesh.update(hexed::criteria::always));
  SECTION("wrong number of BCs") {
    std::vector<std::shared_ptr<hexed::Flow_bc>> bcs;
    for (int i = 0; i < 4; ++i) bcs.emplace_back(new hexed::Copy);
    REQUIRE_THROWS(mesh.add_tree(bcs));
  }
  {
    std::vector<std::shared_ptr<hexed::Flow_bc>> bcs;
    for (int i = 0; i < 6; ++i) bcs.emplace_back(new hexed::Copy);
    mesh.add_tree(bcs);
  }
  SECTION("multiple `add_tree` calls") {
    std::vector<std::shared_ptr<hexed::Flow_bc>> bcs;
    for (int i = 0; i < 6; ++i) bcs.emplace_back(new hexed::Copy);
    REQUIRE_THROWS(mesh.add_tree(bcs));
  }
  REQUIRE(mesh.elements().size() == 1);
  mesh.valid().assert_valid();
  SECTION("refinement") {
    mesh.update();
    REQUIRE(mesh.elements().size() == 8);
    mesh.valid().assert_valid();
    mesh.update([](hexed::Element& elem){
      auto np = elem.nominal_position();
      return np[0] == 0 && np[1] == 0 && np[2] == 0;
    });
    REQUIRE(mesh.elements().size() == 15);
    mesh.valid().assert_valid();
    // refining this element should refine 3 face neighbors and 3 edge neighbors
    mesh.update([](hexed::Element& elem){
      auto np = elem.nominal_position();
      return elem.refinement_level() == 2 && np[0] == 1 && np[1] == 1 && np[2] == 1;
    });
    REQUIRE(mesh.elements().size() == 64);
    mesh.valid().assert_valid();
    mesh.update([](hexed::Element& elem){
      auto np = elem.nominal_position();
      return elem.refinement_level() == 1 && np[0] == 1 && np[1] == 1 && np[2] == 1;
    });
    mesh.update([](hexed::Element& elem){
      auto np = elem.nominal_position();
      return elem.refinement_level() == 2 && np[0] == 2 && np[1] == 2 && np[2] == 2;
    });
    mesh.valid().assert_valid();
    REQUIRE(mesh.elements().size() == 78);
  }
  SECTION("unrefinement") {
    for (int i = 0; i < 3; ++i) mesh.update();
    REQUIRE(mesh.elements().size() == 512);
    mesh.valid().assert_valid();
    auto predicate = [](hexed::Element& elem) {
      auto np = elem.nominal_position();
      int thresh = hexed::math::pow(2, elem.refinement_level())/2;
      return    (np[0] <  thresh && np[1] <  thresh && np[2] <  thresh)
             || (np[0] >= thresh && np[1] >= thresh && np[2] >= thresh);
    };
    mesh.set_unref_locks(hexed::criteria::always);
    // this should do nothing because we just locked unrefinement for all the elements
    mesh.update(hexed::criteria::never, predicate);
    REQUIRE(mesh.elements().size() == 512);
    mesh.set_unref_locks(); // no longer locked
    // now unrefinement should work
    mesh.update(hexed::criteria::never, predicate);
    REQUIRE(mesh.elements().size() == 6*64 + 2*8);
    mesh.valid().assert_valid();
    // this should do nothing because of ref level smoothing
    mesh.update(hexed::criteria::never, predicate);
    REQUIRE(mesh.elements().size() == 6*64 + 2*8);
    mesh.valid().assert_valid();
    // simultaneous refinement and unrefinement
    mesh.update([](hexed::Element& elem){return elem.refinement_level() == 2;},
                [](hexed::Element& elem){return elem.refinement_level() == 3;});
    REQUIRE(mesh.elements().size() == 6*8 + 2*64);
    mesh.valid().assert_valid();
    SECTION("neighbors with different ref levels") {
      hexed::Accessible_mesh mesh1({1, 5, 3, hexed::config::max_row_size}, .7, hexed::laminar);
      mesh1.add_boundary_condition(std::make_shared<hexed::Copy>());
      std::vector<std::shared_ptr<hexed::Flow_bc>> bcs;
      for (int i = 0; i < 6; ++i) bcs.push_back(std::make_shared<hexed::Copy>());
      mesh1.add_tree(bcs);
      mesh1.update();
      mesh1.update();
      mesh1.update([](hexed::Element& elem){return elem.nominal_position()[0] < 2;});
      REQUIRE(mesh1.elements().size() == 4*8 + 4*64);
      mesh1.valid().assert_valid();
      mesh1.update(hexed::criteria::never, hexed::criteria::always);
      REQUIRE(mesh1.elements().size() == 4*1 + 4*8);
      mesh1.valid().assert_valid();
    }
  }
  SECTION("no diagonally-connected elements") {
    mesh.update();
    mesh.update();
    std::vector<hexed::Mat<3, 3>> triangles(2);
    triangles[0] << .7/8, .7/8, .7/8,
                    .7/8, .7/8, .7/8,
                      0.,   0.,   .7;
    triangles[1] << 2.1/8, 2.1/8, 2.1/8,
                    2.1/8, 2.1/8, 2.1/8,
                       0.,    0.,    .7;
    mesh.set_surface(new hexed::Simplex_geom<3>(triangles), std::make_shared<hexed::Copy>(), hexed::Mat<3>{.6, .6, .6});
    // count number of non-extruded elements
    int count = 0;
    auto& elems = mesh.elements();
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      if (elems[i_elem].tree) ++count;
    }
    // number should indicate that the diagonally-connected elements have been deleted
    REQUIRE(count == 48);
  }
}

TEST_CASE("mesh I/O") {
  hexed::printers::info.printers.clear();
  hexed::Mat<3> correct_sum_vertices = hexed::Mat<3>::Zero();
  int correct_n_vertices = 0;
  int correct_n_car_after = 0;
  int correct_n_def_after = 0;
  { // create a mesh and write it to a file
    hexed::Accessible_mesh mesh({1, 4, 2, hexed::config::max_row_size - 1}, .8, hexed::laminar);
    std::vector<std::shared_ptr<hexed::Flow_bc>> bcs;
    for (int i = 0; i < 4; ++i) bcs.emplace_back(new hexed::Copy);
    mesh.add_tree(bcs, hexed::Mat<2>{0.1, 0.2});
    mesh.update();
    mesh.update([](hexed::Element& elem){return elem.nominal_position()[0] != elem.nominal_position()[1];});
    mesh.set_surface(new hexed::Hypersphere(hexed::Mat<2>{.9, 0.2}, 0.1), std::make_shared<hexed::Nonpenetration>());
    REQUIRE(mesh.cartesian().elements().size() == 6);
    REQUIRE(mesh.deformed().elements().size() == 7);
    mesh.write("io_test");
    mesh.visualize("default", "io_test_orig");
    // compute the sum of the vertex coordinates of all elements (counting each vertex once for each element using it)
    // to check vertex position
    auto& elems = mesh.elements();
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      for (int i_vert = 0; i_vert < 4; ++i_vert) {
        correct_sum_vertices += elems[i_elem].shape().vertex(i_vert).point({});
        REQUIRE_THAT(elems[i_elem].shape().vertex(i_vert).point({}), Catch::Matchers::RangeEquals(
                     elems[i_elem].active_shape().vertex(i_vert).point({}), hexed::math::Approx_equal(0., 1e-8)));
        ++correct_n_vertices;
      }
    }
    // refine the mesh again and count the number of Cartesian and deformed elements
    // to make sure the recreated mesh behaves the same way
    mesh.update([](hexed::Element& elem){return elem.nominal_position()[0] > 2;});
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      if (elems[i_elem].get_is_deformed()) ++correct_n_def_after;
      else                                 ++correct_n_car_after;
    }
  }
  { // read the above mesh from the file and check that it's the same
    std::vector<std::shared_ptr<hexed::Flow_bc>> extr_bcs;
    for (int i = 0; i < 4; ++i) extr_bcs.emplace_back(new hexed::Copy);
    hexed::Accessible_mesh mesh("io_test", extr_bcs, hexed::laminar,
                                new hexed::Hypersphere(hexed::Mat<2>{.9, 0.2}, 0.1),
                                std::make_shared<hexed::Nonpenetration>());
    mesh.visualize("default", "io_test_reconstructed");
    REQUIRE(mesh.root_size() == Catch::Approx(0.8));
    REQUIRE(mesh.cartesian().elements().size() == 6);
    REQUIRE(mesh.deformed().elements().size() == 7);
    auto& elems = mesh.elements();
    int rl1 = 0;
    int rl2 = 0;
    hexed::Mat<3> sum_vertices = hexed::Mat<3>::Zero();
    int n_vertices = 0;
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      rl1 += elems[i_elem].refinement_level() == 1;
      rl2 += elems[i_elem].refinement_level() == 2;
      REQUIRE(elems[i_elem].tree);
      for (int i_vert = 0; i_vert < 4; ++i_vert) {
        sum_vertices += elems[i_elem].shape().vertex(i_vert).point({});
        REQUIRE_THAT(elems[i_elem].shape().vertex(i_vert).point({}), Catch::Matchers::RangeEquals(
                     elems[i_elem].active_shape().vertex(i_vert).point({}), hexed::math::Approx_equal(0., 1e-8)));
        ++n_vertices;
      }
    }
    REQUIRE(rl1 == 2);
    REQUIRE(rl2 == 11);
    REQUIRE(n_vertices == correct_n_vertices);
    REQUIRE(sum_vertices(0) == Catch::Approx(correct_sum_vertices(0)));
    REQUIRE(sum_vertices(1) == Catch::Approx(correct_sum_vertices(1)));
    mesh.valid().assert_valid();
    // refine the mesh and check that it's the same as refining the original mesh
    mesh.update([](hexed::Element& elem){return elem.nominal_position()[0] > 2;});
    int n_car_after = 0;
    int n_def_after = 0;
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      if (elems[i_elem].get_is_deformed()) ++n_def_after;
      else                                 ++n_car_after;
    }
    mesh.valid().assert_valid();
    REQUIRE(n_car_after == correct_n_car_after);
    REQUIRE(n_def_after == correct_n_def_after);
  }
}

TEST_CASE("masking") {
  hexed::Accessible_mesh mesh({2, 4, 2, 2}, 1., hexed::laminar);
  hexed::Gauss_legendre basis(2);
  std::vector<std::shared_ptr<hexed::Flow_bc>> bcs;
  for (int i = 0; i < 4; ++i) bcs.emplace_back(new hexed::Copy);
  mesh.add_tree(bcs);
  mesh.update();
  mesh.update([](hexed::Element& elem){return elem.nominal_position()[0] == elem.nominal_position()[1];});
  SECTION("initial mask") {
    mesh.reset_masks();
    hexed::Accessible_mesh::Masked_mesh(mesh, basis);
    hexed::Accessible_mesh::Masked_mesh masked(mesh, basis);
    auto& elems = mesh.elements();
    REQUIRE(masked.kernel_mesh.n_dim == 2);
    REQUIRE(masked.kernel_mesh.row_size == 2);
    int n_masked = 0;
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) n_masked += elems[i_elem].mask();
    REQUIRE(n_masked == 10);
    CHECK(masked.kernel_mesh.elems.size() == 10);
    CHECK(masked.kernel_mesh.car_connections.size() == 28);
    CHECK(masked.kernel_mesh.def_connections.size() == 0);
    CHECK(masked.kernel_mesh.face_refinements.size() == 4);
    CHECK(masked.bound_cons.size()  == 12);
  }
  SECTION("custom mask") {
    mesh.reset_masks();
    hexed::Accessible_mesh::Masked_mesh(mesh, basis);
    hexed::Accessible_mesh::Masked_mesh masked(mesh, basis, [](hexed::Element& elem){
      return elem.shape().vertex(3).point({})[1] < .501;
    });
    auto& elems = mesh.elements();
    int n_masked = 0;
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      n_masked += elems[i_elem].mask();
    }
    REQUIRE(n_masked == 5);
    REQUIRE(masked.kernel_mesh.elems.size() == 5);
    REQUIRE(masked.kernel_mesh.car_elems.size() == 5);
    REQUIRE(masked.kernel_mesh.def_elems.size() == 0);
    REQUIRE(masked.kernel_mesh.car_connections.size() == 16);
    REQUIRE(masked.kernel_mesh.def_connections.size() == 0);
    REQUIRE(masked.kernel_mesh.face_refinements.size() == 3);
    REQUIRE(masked.bound_cons.size() == 6);
  }
}
