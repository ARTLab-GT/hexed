#include <catch2/catch.hpp>
#include <Vis_data.hpp>
#include <Deformed_element.hpp>
#include <Domain_func.hpp>
#include <Gauss_legendre.hpp>
#include <config.hpp>
#include <Spacetime_func.hpp>

class Rad_sq : public hexed::Spacetime_func
{
  Eigen::Vector3d center;
  public:
  Rad_sq(Eigen::Vector3d c) : center{c} {}
  virtual int n_var(int n_dim) const {return 1;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    Eigen::Map<Eigen::VectorXd> p(pos.data(), pos.size());
    Eigen::Vector3d q = p - center;
    return {(q.transpose()*q)[0]};
  }
};

TEST_CASE("Vis_data")
{
  // create some arbitrarily-shaped elements
  hexed::Deformed_element elem3({2, 5, 3, hexed::config::max_row_size});
  hexed::Deformed_element elem2({2, 4, 2, hexed::config::max_row_size});
  hexed::Gauss_legendre basis(hexed::config::max_row_size);
  elem3.vertex(7).pos = {.9, .9, .9};
  elem2.vertex(2).pos = {.8, .3, 0.};
  hexed::Vis_data vis3(elem3, hexed::Position_func(), basis);
  hexed::Vis_data vis2(elem2, hexed::Position_func(), basis);
  hexed::Vis_data vis_const(elem2, hexed::Constant_func({0.3, 0.4, 0.5, 0.6}), basis);

  SECTION("edges")
  {
    SECTION("2D")
    {
      auto edges = vis2.edges();
      REQUIRE(edges.size() == 21*2*4);
      REQUIRE(edges( 0 + 21*2*0) == Approx(0.).scale(1.));
      REQUIRE(edges( 1 + 21*2*2) == Approx(0.).scale(1.));
      REQUIRE(edges(22 + 21*2*2) == Approx(.05).scale(1.));
      REQUIRE(edges( 1 + 21*2*0) == Approx(.8*.05).scale(1.));
      REQUIRE(edges(22 + 21*2*0) == Approx(.3*.05).scale(1.));
    }
    SECTION("3D")
    {
      auto edges = vis3.edges(6);
      REQUIRE(edges.size() == 6*3*12);
      REQUIRE(edges( 0 + 6*3*0) == Approx(0.).scale(1.));
      REQUIRE(edges( 1 + 6*3*0) == Approx(0.2).scale(1.));
      REQUIRE(edges( 7 + 6*3*0) == Approx(0.).scale(1.));
      REQUIRE(edges(13 + 6*3*0) == Approx(0.).scale(1.));
      REQUIRE(edges( 1 + 6*3*1) == Approx(0.2).scale(1.));
      REQUIRE(edges( 7 + 6*3*1) == Approx(0.).scale(1.));
      REQUIRE(edges(13 + 6*3*1) == Approx(1.).scale(1.));
      REQUIRE(edges( 1 + 6*3*4) == Approx(0.).scale(1.));
      REQUIRE(edges( 7 + 6*3*4) == Approx(0.2).scale(1.));
      REQUIRE(edges(13 + 6*3*4) == Approx(0.).scale(1.));
      REQUIRE(edges( 2 + 6*3*6) == Approx(1.).scale(1.));
      REQUIRE(edges( 8 + 6*3*6) == Approx(0.4).scale(1.));
      REQUIRE(edges(14 + 6*3*6) == Approx(0.).scale(1.));
      REQUIRE(edges( 5 + 6*3*11) == Approx(.9).scale(1.));
      REQUIRE(edges(11 + 6*3*11) == Approx(.9).scale(1.));
      REQUIRE(edges(17 + 6*3*11) == Approx(.9).scale(1.));
    }
    SECTION("constant")
    {
      REQUIRE(vis_const.edges().size() == 21*4*4);
    }
  }

  SECTION("interior")
  {
    SECTION("2D")
    {
      auto interior = vis2.interior();
      REQUIRE(interior.size() == 21*21*2);
      REQUIRE(interior(0) == Approx(0.).scale(1.));
      REQUIRE(interior(         1) == Approx(0.).scale(1.));
      REQUIRE(interior(21*21 +  1) == Approx(.05).scale(1.));
      REQUIRE(interior(        21) == Approx(.05*.8).scale(1.));
      REQUIRE(interior(21*21 + 21) == Approx(.05*.3).scale(1.));
      REQUIRE(interior(21*21*1 - 1) == Approx(1.).scale(1.));
      REQUIRE(interior(21*21*2 - 1) == Approx(1.).scale(1.));
    }
    SECTION("3D")
    {
      auto interior = vis3.interior();
      REQUIRE(interior.size() == 21*21*21*3);
      REQUIRE(interior(21*21*21*0) == Approx(0.).scale(1.));
      REQUIRE(interior(21*21*21*1) == Approx(0.).scale(1.));
      REQUIRE(interior(21*21*21*2) == Approx(0.).scale(1.));
      REQUIRE(interior(21*21*21*1 - 1) == Approx(.9).scale(1.));
      REQUIRE(interior(21*21*21*2 - 1) == Approx(.9).scale(1.));
      REQUIRE(interior(21*21*21*3 - 1) == Approx(.9).scale(1.));
    }
    SECTION("constant")
    {
      auto interior = vis_const.interior();
      REQUIRE(interior.size() == 21*21*4);
      REQUIRE(interior(21*21*0) == Approx(.3).scale(1.));
      REQUIRE(interior(21*21*1) == Approx(.4).scale(1.));
      REQUIRE(interior(21*21*2) == Approx(.5).scale(1.));
      REQUIRE(interior(21*21*3) == Approx(.6).scale(1.));
    }
  }

  SECTION("face")
  {
    auto face = vis2.face(1, 0);
    REQUIRE(face.size() == 21*2);
    REQUIRE(face(0*21 + 0) == Approx(0.).scale(1.));
    REQUIRE(face(1*21 + 0) == Approx(0.).scale(1.));
    REQUIRE(face(0*21 + 2) == Approx(2*.05*.8).scale(1.));
    REQUIRE(face(1*21 + 2) == Approx(2*.05*.3).scale(1.));
  }

  SECTION("sample")
  {
    auto sample = vis3.sample(Eigen::Matrix<double, 3, 2>{{.5, 1.}, {.5, 0.}, {.5, 0.}});
    REQUIRE(sample.rows() == 3);
    REQUIRE(sample.cols() == 2);
    REQUIRE((sample(Eigen::all, 0) - Eigen::VectorXd::Constant(3, .5 - .1/8.)).norm() == Approx(0.).scale(1.));
    REQUIRE((sample(Eigen::all, 1) - Eigen::Vector3d{1., 0., 0.}).norm() == Approx(0.).scale(1.));
  }

  SECTION("contour")
  {
    SECTION("corner")
    {
      hexed::Element elem({2, 5, 3, hexed::config::max_row_size});
      Rad_sq func({0, 1, 0});
      hexed::Vis_data vis(elem, func, basis);
      auto con = vis.compute_contour(.04, 2); // exact contour is quarter-sphere with radius .2 centered at (0, 1, 0)
      REQUIRE(con.vert_ref_coords.rows() == 7);
      REQUIRE(con.vert_ref_coords.cols() == 3);
      REQUIRE(con.normals.rows() == 7);
      REQUIRE(con.normals.cols() == 3);
      // don't test number of quads because there are some degenerate ones
      REQUIRE(con.elem_vert_inds.cols() == 4);
      // test vertices are approximately on the correct contour surface
      for (int i_vert = 0; i_vert < con.vert_ref_coords.rows(); ++i_vert) {
        std::vector<double> pos(3);
        for (int i_dim = 0; i_dim < 3; ++i_dim) pos[i_dim] = con.vert_ref_coords(i_vert, i_dim);
        REQUIRE(func(pos, 0.)[0] == Approx(.04).margin(1e-6));
        Eigen::MatrixXd rad = con.vert_ref_coords(i_vert, Eigen::all) - Eigen::Matrix<double, 1, 3>{0, 1, 0};
        REQUIRE((con.normals(i_vert, Eigen::all) - rad.normalized()).norm() == Approx(0.).scale(1.));
      }
    }
    SECTION("center")
    {
      hexed::Element elem({2, 5, 3, hexed::config::max_row_size});
      hexed::Vis_data vis(elem, Rad_sq({.5, .5, .5}), basis);
      auto con = vis.compute_contour(.04, 2); // exact contour is sphere with radius .2 centered at (.5, .5, .5)
      REQUIRE(con.vert_ref_coords.rows() == 26);
      REQUIRE(con.vert_ref_coords.cols() == 3);
      REQUIRE(con.normals.rows() == 26);
      REQUIRE(con.normals.cols() == 3);
      REQUIRE(con.elem_vert_inds.rows() == 24);
      REQUIRE(con.elem_vert_inds.cols() == 4);
    }
    SECTION("warped")
    {
      hexed::Deformed_element elem({2, 4, 2, hexed::config::max_row_size});
      for (int i_vert = 1; i_vert < 4; i_vert += 2) {
        elem.vertex(i_vert).pos[0] += 0.5;
      }
      elem.set_jacobian(basis);
      hexed::Vis_data vis(elem, hexed::Linear(Eigen::Vector2d{1., 0.}), basis);
      auto con = vis.compute_contour(.5, 3);
      for (int i_vert = 0; i_vert < con.vert_ref_coords.rows(); ++i_vert) {
        REQUIRE((con.normals(i_vert, Eigen::all) - Eigen::Matrix<double, 1, 2>{1., 0.}).norm() == Approx(0.).scale(1.));
      }
    }
  }
}
