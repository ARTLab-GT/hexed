#include <hexed/Simplex_geom.hpp>
#include <hexed/Tecplot_file.hpp>
#include <hexed/Xdmf_wrapper.hpp>
#include <hexed/Spacetime_func.hpp>

namespace hexed {

#if HEXED_OBSESSIVE_TIMING
Stopwatch_tree Simplex_geom_nd::stopwatch("", {
  {"nearest_point", Stopwatch_tree("projection")},
  {"intersections", Stopwatch_tree("intersection")}
});
#endif

template<>
void Simplex_geom<2>::merge(Nearest_point<2>& nearest, Mat<2, 2> sim, Mat<2> point) {
  nearest.merge(math::proj_to_segment({sim(all, 0), sim(all, 1)}, point));
}

template<>
void Simplex_geom<3>::merge(Nearest_point<3>& nearest, Mat<3, 3> sim, Mat<3> point) {
  // try projecting the point to the plane of the triangle
  Mat<3, 2> lhs;
  lhs(all, 0) = sim(all, 1) - sim(all, 0);
  lhs(all, 1) = sim(all, 2) - sim(all, 0);
  Mat<2> lstsq = (lhs.transpose()*lhs).inverse()*(lhs.transpose()*(point - sim(all, 0)));
  if (lstsq(0) >= 0 && lstsq(1) >= 0 && lstsq.sum() <= 1) {
    // if the projected point is inside the triangle, evaluate it as the potential nearest point
    nearest.merge(sim(all, 0) + lhs*lstsq);
  } else {
    // if the projected point is outside the triangle,
    // fall back to finding the nearest point on all the edges of the triangle
    for (int i_edge = 0; i_edge < 3; ++i_edge) {
      nearest.merge(math::proj_to_segment({sim(all, i_edge), sim(all, (i_edge + 1)%3)}, point));
    }
  }
}

//! \cond
template<>
void Simplex_geom<3>::visualize(std::string format, std::string fname) {
  Array<int> triangles({int(_simplices.size()), 3});
  Array<double> pos({3, 3*int(_simplices.size())});
  Array<double> vars({0, 3*int(_simplices.size())});
  int i_vert = 0;
  for (Mat<3, 3> sim : _simplices) {
    for (int elem_vert = 0; elem_vert < 3; ++elem_vert, ++i_vert) {
      triangles[i_vert] = i_vert;
      for (int i_dim = 0; i_dim < 3; ++i_dim) pos(i_dim)[i_vert] = sim(i_dim, elem_vert);
    }
  }
  Visualizer::create(format, 3, 2, fname, {}, 0., Visualizer::simplex)->write_unstruct(triangles, pos, vars);
  auto edge_vis = Visualizer::create(format, 3, 1, fname + "_edges", {}, 0., Visualizer::block);
  for (auto& edge : _geom_edges) {
    Array<double> points = edge.points();
    Array<double> transposed({3, edge.n_points()});
    for (int i_point = 0; i_point < edge.n_points(); ++i_point) {
      for (int i_dim = 0; i_dim < 3; ++i_dim) transposed(i_dim)[i_point] = points(i_point)[i_dim];
    }
    edge_vis->write_block(transposed, Array<double>({0, edge.n_points()}));
  }
}
//! \endcond

std::vector<Mat<2, 2>> segments(const Mat<dyn, dyn>& points) {
  std::vector<Mat<2, 2>> sims;
  for (unsigned i = 0; i < points.cols() - 1; ++i) {
    Mat<2, 2> sim;
    sim(all, 0) = points(all, i);
    sim(all, 1) = points(all, i + 1);
    sims.push_back(sim);
  }
  return sims;
}

}
