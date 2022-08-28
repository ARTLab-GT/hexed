#include <otter_vis.hpp>
#if CARTDG_USE_OTTER
#include <otter/points.hpp>
#include <Vis_data.hpp>
#include <Spacetime_func.hpp>

namespace cartdg::otter_vis
{

void add_edges(otter::plot& plt, Element& elem, const Basis& bas, int n_sample)
{
  const int n_dim = elem.storage_params().n_dim;
  Vis_data vis(elem, Position_func(), bas);
  Eigen::MatrixXd edges = vis.edges(n_sample);
  edges.resize(n_sample, edges.size()/n_sample);
  for (int i_edge = 0; i_edge < edges.cols()/n_dim; ++i_edge) {
    plt.add(otter::curve(edges(Eigen::all, Eigen::seqN(i_edge*n_dim, n_dim))));
  }
}

void add_contour(otter::plot& plt, Element& elem, const Basis& bas, const Qpoint_func& func,
                 double value, int n_div, color_spec spec, double time, int i_var)
{
  Vis_data vis(elem, func, bas, time);
  auto con = vis.compute_contour(value, n_div, i_var);
  Vis_data vis_pos(elem, Position_func(), bas, time);
  Vis_data vis_col(elem, spec.color_by, bas, time);
  if (elem.storage_params().n_dim == 3) {
    otter::surface_data data;
    for (int i_vert = 0; i_vert < con.vert_ref_coords.rows(); ++i_vert) {
      otter::surface_data::vertex vert;
      vert.pos = vis_pos.sample(con.vert_ref_coords(i_vert, Eigen::all).transpose());
      vert.normal = con.normals(i_vert, Eigen::all).transpose();
      Eigen::MatrixXd color = spec.map((vis_col.sample(vert.pos)(0) - spec.bounds[0])/(spec.bounds[1] - spec.bounds[0]));
      vert.rgba(Eigen::seq(0, color.size())) = color.transpose();
      data.verts.push_back(vert);
    }
    for (int i_quad = 0; i_quad < con.elem_vert_inds.rows(); ++i_quad) {
      data.add_quad(con.elem_vert_inds(i_quad, Eigen::all).transpose());
    }
    otter::surface surface(data);
    otter::surface flipped(surface);
    flipped.flip_orientation();
    surface.append(flipped);
    plt.add(surface);
  }
}

}
#endif
