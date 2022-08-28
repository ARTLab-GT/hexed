#include <Vis_data.hpp>
#include <math.hpp>
#include <otter/plot.hpp>
#include <otter/colors.hpp>
#include <iostream>

namespace cartdg
{

Eigen::MatrixXd Vis_data::sample_qpoint_data(Eigen::VectorXd qpoint_data, Eigen::MatrixXd ref_coords)
{
  int nv = qpoint_data.size()/n_qpoint;
  const int n_sample = ref_coords.cols();
  Eigen::MatrixXd result(nv, n_sample);
  for (int i_sample = 0; i_sample < n_sample; ++i_sample) {
    // ith row is the interpolation matrix along the ith dimension
    auto interp = bas.interpolate(ref_coords(Eigen::all, i_sample));
    for (int i_var = 0; i_var < nv; ++i_var) {
      // start with all the data for this variable
      Eigen::VectorXd var = qpoint_data(Eigen::seqN(i_var*n_qpoint, n_qpoint));
      // interpolate one dimension at a time
      for (int i_dim = n_dim - 1; i_dim >= 0; --i_dim) {
        var = custom_math::dimension_matvec(interp(i_dim, Eigen::all), var, i_dim);
      }
      // ...until all you have left is a vector with one element
      result(i_var, i_sample) = var(0);
    }
  }
  return result;
}

Vis_data::Vis_data(Element& elem, const Qpoint_func& func, const Basis& basis, double time) :
  n_dim{elem.storage_params().n_dim},
  n_edge{custom_math::pow(2, n_dim - 1)*n_dim},
  row_size{elem.storage_params().row_size},
  n_qpoint{elem.storage_params().n_qpoint()},
  n_var{func.n_var(n_dim)},
  bas{basis},
  vars(n_qpoint*n_var)
{
  // fetch data at quadrature points
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    auto v = func(elem, bas, i_qpoint, time);
    for (int i_var = 0; i_var < n_var; ++i_var) vars(n_qpoint*i_var + i_qpoint) = v[i_var];
  }
}

Eigen::VectorXd Vis_data::edges(int n_sample)
{
  Eigen::MatrixXd interp {bas.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  const int nfqpoint = n_qpoint/row_size;
  Eigen::MatrixXd result(n_var*n_sample, n_edge);
  // interpolate qpoint varsition to edges
  Eigen::MatrixXd boundary {bas.boundary()};
  // interpolate all edges which point in direction `i_dim`
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    const int stride {custom_math::pow(row_size, n_dim - 1 - i_dim)};
    const int n_outer {n_qpoint/stride/row_size};
    // extract the `i_var`th variable
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      Eigen::MatrixXd edge_qpoints {row_size, n_edge/n_dim}; // quadrature points interpolated to the edge
      for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
        Eigen::VectorXd qpoint_slab {nfqpoint};
        for (int i_outer = 0; i_outer < n_outer; ++i_outer) {
          for (int i_inner = 0; i_inner < stride; ++i_inner) {
            qpoint_slab[i_outer*stride + i_inner] = vars(i_var*n_qpoint + i_qpoint*stride + i_outer*stride*row_size + i_inner);
          }
        }
        // interpolate edge quadrature points to edge visualization points
        edge_qpoints.row(i_qpoint) = custom_math::hypercube_matvec(boundary, qpoint_slab);
      }
      result(Eigen::seqN(i_var*n_sample, n_sample),
             Eigen::seqN(i_dim*n_edge/n_dim, n_edge/n_dim)) = interp*edge_qpoints; // note: rhs auto resized to fit lhs
    }
  }
  result.resize(result.size(), 1);
  return result;
}

Eigen::VectorXd Vis_data::interior(int n_sample)
{
  Eigen::MatrixXd interp {bas.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  const int n_block = custom_math::pow(n_sample, n_dim);
  Eigen::VectorXd result(n_block*n_var);
  for (int i_var = 0; i_var < n_var; ++i_var) {
    result(Eigen::seqN(i_var*n_block, n_block)) = custom_math::hypercube_matvec(interp, vars(Eigen::seqN(i_var*n_qpoint, n_qpoint)));
  }
  return result;
}

Eigen::VectorXd Vis_data::face(int i_dim, bool is_positive, int n_sample)
{
  Eigen::MatrixXd interp {bas.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  Eigen::MatrixXd bound = bas.boundary()(is_positive, Eigen::all);
  const int n_block = custom_math::pow(n_sample, n_dim - 1);
  Eigen::VectorXd result(n_block*n_var);
  for (int i_var = 0; i_var < n_var; ++i_var) {
    // interpolate quadrature points to face quadrature points
    auto var = vars(Eigen::seqN(i_var*n_qpoint, n_qpoint));
    Eigen::VectorXd face_qpoints = custom_math::dimension_matvec(bound, var, i_dim);
    // interpolate face quadrature points to uniformly spaced
    result(Eigen::seqN(i_var*n_block, n_block)) = custom_math::hypercube_matvec(interp, face_qpoints);
  }
  return result;
}

Eigen::MatrixXd Vis_data::sample(Eigen::MatrixXd ref_coords)
{
  return sample_qpoint_data(vars, ref_coords);
}

Vis_data::Contour Vis_data::compute_contour(double value, int n_div, int i_var)
{
  otter::plot plt;
  Contour con;
  auto sample = interior(n_div + 1);
  const int n_sample = sample.size();
  const int n_block = custom_math::pow(2*n_div + 1, n_dim);
  const int n_corner = custom_math::pow(2, n_dim - 1);
  Eigen::VectorXi i_contour = Eigen::VectorXi::Constant(n_block, -1);
  std::vector<int> i_block;
  std::vector<int> faces;
  std::vector<int> strides_sample;
  std::vector<int> strides_block;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    strides_sample.push_back(custom_math::pow(n_div + 1, n_dim - 1 - i_dim));
    strides_block.push_back(custom_math::pow(2*n_div + 1, n_dim - 1 - i_dim));
  }
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    int stride = custom_math::pow(n_div + 1, n_dim - i_dim - 1);
    for (int i_outer = 0; i_outer < n_sample/stride/(n_div + 1); ++i_outer) {
      for (int i_inner = 0; i_inner < stride; ++i_inner) {
        for (int i_row = 0; i_row < n_div; ++i_row) { // don't do the last point
          int sample0 = (i_outer*(n_div + 1) + i_row)*stride + i_inner;
          int sample1 = sample0 + stride;
          Eigen::VectorXd coords(n_dim);
          for (int j_dim = 0; j_dim < n_dim; ++j_dim) coords(j_dim) = ((sample0/strides_sample[j_dim])%(n_div + 1))/double(n_div);
          plt.add(otter::points(coords.transpose(), otter::colors::tableau[1]));
          if ((sample[sample0] > value) != (sample[sample1] > value)) {
            int block = strides_block[i_dim];
            bool boundary [3][2];
            for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
              int row = (sample0/strides_sample[j_dim])%(n_div + 1);
              boundary[j_dim][0] = row == 0;
              boundary[j_dim][1] = row == n_div;
              block += row*2*strides_block[j_dim];
            }
            int n_vert = custom_math::pow(3, n_dim - 1);
            std::vector<int> verts(n_vert);
            for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
              int vert = block;
              for (int j_dim = 0, k_dim = 0; j_dim < n_dim; ++j_dim) if (j_dim != i_dim) {
                int row_incr = ((i_vert/custom_math::pow(3, n_dim - 2 - k_dim))%3 - 1);
                vert += (!boundary[j_dim][row_incr > 0])*((i_vert/custom_math::pow(3, n_dim - 2 - k_dim))%3 - 1)*strides_block[j_dim];
                ++k_dim;
              }
              int ib;
              for (int j_dim = 0; j_dim < n_dim; ++j_dim) coords(j_dim) = ((vert/strides_block[j_dim])%(2*n_div + 1))/(2.*n_div);
              plt.add(otter::points(coords.transpose(), otter::colors::tableau[2]));
              if (i_contour(vert) < 0) {
                ib = i_block.size();
                i_block.push_back(vert);
                i_contour(vert) = ib;
              } else ib = i_contour(vert);
              verts[i_vert] = ib;
            }
            for (int i_elem = 0; i_elem < n_corner; ++i_elem) {
              for (int i_corner = 0; i_corner < n_corner; ++i_corner) {
                int i_vert = 0;
                for (int j_dim = 0; j_dim < n_dim - 1; ++j_dim) {
                  int stride = custom_math::pow(2, n_dim - 2 - j_dim);
                  i_vert += ((i_elem/stride)%2 + (i_corner/stride)%2)*custom_math::pow(3, n_dim - 2 - j_dim);
                }
                faces.push_back(verts[i_vert]);
              }
            }
          }
        }
      }
    }
  }
  con.vert_ref_coords.resize(i_block.size(), n_dim);
  for (unsigned i = 0; i < i_block.size(); ++i) {
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      con.vert_ref_coords(i, i_dim) = ((i_block[i]/strides_block[i_dim])%(2*n_div + 1))/(2.*n_div);
    }
  }
  con.normals.resize(i_block.size(), n_dim);
  con.normals.setZero();
  con.elem_vert_inds.resize(faces.size()/n_corner, n_corner);
  for (int i_elem = 0; i_elem < int(faces.size())/n_corner; ++i_elem) {
    for (int i_corner = 0; i_corner < n_corner; ++i_corner) {
      con.elem_vert_inds(i_elem, i_corner) = faces[i_elem*n_corner + i_corner];
    }
  }
  Eigen::MatrixXd verts(8, 3);
  for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) for (int k = 0; k < 2; ++k) {
      verts(k + 2*(j + 2*i), Eigen::all) << i, j, k;
  }
  plt.add(otter::points(verts, otter::colors::tableau[0]));
  plt.show();
  return con;
}

}
