#include <hexed/Vis_data.hpp>

namespace hexed {

Array<double> check(Array<double>&& data) {
  HEXED_ASSERT(data.order() >= 2, "Order of `data` is too small.");
  auto shape = data.shape();
  for (int i_dim = 2; i_dim < data.order(); ++i_dim) {
    HEXED_ASSERT(shape[i_dim] == shape[1], "All dimensions except the first must be the same.");
  }
  return data;
}

Vis_data::Vis_data(Array<double> data, const Basis& basis)
: _data{check(std::move(data))}
, _n_dim{_data.order() - 1}
, _n_edge{math::pow(2, _n_dim - 1)*_n_dim}
, _row_size{_data.shape()[1]}
, _n_qpoint{_data(0).size()}
, _n_var{_data.shape()[0]}
, _basis{basis}
{}

Array<double> Vis_data::sample(Array<double> coords) const {
  HEXED_ASSERT(coords.order() >= 2, "`coords` must be order-2.");
  HEXED_ASSERT(coords.shape()[0] == _n_dim, "`coords` must have one row for each reference coordinate.");
  Array<double> s({_n_var, coords(0).size()});
  for (int i_var = 0; i_var < _n_var; ++i_var) s(i_var) = _sample(_data(i_var), coords());
  return s;
}

Array<double> Vis_data::interior(int n_sample) const {
  Array<double> result(hypercubes(_n_var, _n_dim, n_sample));
  Mat<dyn, dyn> interp = _basis.interpolate(Mat<>::LinSpaced(n_sample, 0., 1.));
  for (int i_var = 0; i_var < _n_var; ++i_var) {
    result(i_var).vector() = math::hypercube_matvec(interp, _data(i_var).vector());
  }
  return result;
}

Array<double> Vis_data::_sample(Array<double> data, Array<double> coords) const {
  Int n_sample = coords(0).size();
  Array<double> s({n_sample});
  for (int i_sample = 0; i_sample < n_sample; ++i_sample) {
    // ith row is the interpolation matrix along the ith dimension
    auto interp = _basis.interpolate(coords.column(i_sample).vector());
    Eigen::VectorXd var = data.vector();
    // interpolate one dimension at a time
    for (int i_dim = _n_dim - 1; i_dim >= 0; --i_dim) {
      var = math::dimension_matvec(interp(i_dim, Eigen::all), var, i_dim);
    }
    // ...until all you have left is a vector with one element
    s[i_sample] = var(0);
  }
  return s;
}

#if 0
Eigen::MatrixXd Vis_data::sample_qpoint_data(Eigen::VectorXd qpoint_data, Eigen::MatrixXd ref_coords) {
  int nv = qpoint_data.size()/n_qpoint;
  const int n_sample = ref_coords.rows();
  Eigen::MatrixXd result(n_sample, nv);
  for (int i_sample = 0; i_sample < n_sample; ++i_sample) {
    // ith row is the interpolation matrix along the ith dimension
    auto interp = bas.interpolate(ref_coords(i_sample, all));
    for (int i_var = 0; i_var < nv; ++i_var) {
      // start with all the data for this variable
      Eigen::VectorXd var = qpoint_data(Eigen::seqN(i_var*n_qpoint, n_qpoint));
      // interpolate one dimension at a time
      for (int i_dim = n_dim - 1; i_dim >= 0; --i_dim) {
        var = math::dimension_matvec(interp(i_dim, Eigen::all), var, i_dim);
      }
      // ...until all you have left is a vector with one element
      result(i_sample, i_var) = var(0);
    }
  }
  return result;
}

Eigen::VectorXd Vis_data::edges(int n_sample) {
  Eigen::MatrixXd interp {bas.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  const int nfqpoint = n_qpoint/row_size;
  Eigen::MatrixXd result(n_var*n_sample, n_edge);
  // interpolate qpoint varsition to edges
  Eigen::MatrixXd boundary {bas.boundary()};
  // interpolate all edges which point in direction `i_dim`
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    const int stride {math::pow(row_size, n_dim - 1 - i_dim)};
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
        edge_qpoints.row(i_qpoint) = math::hypercube_matvec(boundary, qpoint_slab);
      }
      result(Eigen::seqN(i_var*n_sample, n_sample),
             Eigen::seqN(i_dim*n_edge/n_dim, n_edge/n_dim)) = interp*edge_qpoints; // note: rhs auto resized to fit lhs
    }
  }
  result.resize(result.size(), 1);
  return result;
}

Eigen::VectorXd Vis_data::interior(int n_sample) {
  Eigen::MatrixXd interp {bas.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  const int n_block = math::pow(n_sample, n_dim);
  Eigen::VectorXd result(n_block*n_var);
  for (int i_var = 0; i_var < n_var; ++i_var) {
    result(Eigen::seqN(i_var*n_block, n_block)) = math::hypercube_matvec(interp, vars(Eigen::seqN(i_var*n_qpoint, n_qpoint)));
  }
  return result;
}

Eigen::VectorXd Vis_data::face(int i_dim, bool is_positive, int n_sample) {
  Eigen::MatrixXd interp {bas.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  Eigen::MatrixXd bound = bas.boundary()(is_positive, Eigen::all);
  const int n_block = math::pow(n_sample, n_dim - 1);
  Eigen::VectorXd result(n_block*n_var);
  for (int i_var = 0; i_var < n_var; ++i_var) {
    // interpolate quadrature points to face quadrature points
    auto var = vars(Eigen::seqN(i_var*n_qpoint, n_qpoint));
    Eigen::VectorXd face_qpoints = math::dimension_matvec(bound, var, i_dim);
    // interpolate face quadrature points to uniformly spaced
    result(Eigen::seqN(i_var*n_block, n_block)) = math::hypercube_matvec(interp, face_qpoints);
  }
  return result;
}

Eigen::MatrixXd Vis_data::sample(Eigen::MatrixXd ref_coords) {
  return sample_qpoint_data(vars, ref_coords);
}

Vis_data::Contour Vis_data::compute_contour(double value, int n_div, int n_newton, double tol) {
  Contour con;
  // sample points used for identifying the contour vertices
  const int n_sample = math::pow(n_div + 1, n_dim);
  Mat<> sample = interior(n_div + 1)(Eigen::seqN((n_var - 1)*n_sample, n_sample));
  // if the candidate vertices that could be in the contour were selected from a
  // uniformly spaced block, how many points would this block have?
  const int n_block = math::pow(2*n_div + 1, n_dim);
  // number of corners of a contour element
  const int n_corner = math::pow(2, n_dim - 1);
  // list of vertices on the contour, by their index in the hypothetical candidate block
  std::vector<int> i_block;
  // for each vertex in the candidate block, what is its index in `i_block`?
  // e.g. the first vertex identified to be on the contour will be 0, the next will be 1, etc.
  // -1 indicates not on the contour (true for most of the vertices)
  Eigen::VectorXi i_contour = Eigen::VectorXi::Constant(n_block, -1);
  std::vector<Eigen::VectorXd> directions; // line search direction for projecting vertices onto contour surface
  std::vector<int> faces; // layout [i_element][i_corner]
  std::vector<int> strides_sample; // strides in the sample block
  std::vector<int> strides_block; // strides in the candidate block
  // compute strides
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    strides_sample.push_back(math::pow(n_div + 1, n_dim - 1 - i_dim));
    strides_block.push_back(math::pow(2*n_div + 1, n_dim - 1 - i_dim));
  }
  // choose vertices to be on contour
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) { // identify all vertices which are on a face in direction +-`i_dim`
    int stride = math::pow(n_div + 1, n_dim - i_dim - 1);
    for (int i_outer = 0; i_outer < n_sample/stride/(n_div + 1); ++i_outer) {
      for (int i_inner = 0; i_inner < stride; ++i_inner) {
        for (int i_row = 0; i_row < n_div; ++i_row) { // don't do the last point
          // the goal is now to decide whether the candidate vertex between `sample0` and `sample1`
          // should be on the contour
          int sample0 = (i_outer*(n_div + 1) + i_row)*stride + i_inner;
          int sample1 = sample0 + stride;
          // if this is true, then the correct contour surface is between `sample0` and `sample1`
          if (((sample[sample0] > value) != (sample[sample1] > value))
              && (std::abs(sample[sample0] - sample[sample1]) > tol)) { // comparison with `tol` avoids spurious contours on constant data
            int block = strides_block[i_dim]; // index of candidate vertex (calculation isn't done yet)
            bool boundary [3][2]; // in a given direction, are we on the boundary? layout: [j_dim][face_positive]
            // evaluate `boundary` and finish calculating `block`
            for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
              int row = (sample0/strides_sample[j_dim])%(n_div + 1);
              boundary[j_dim][0] = row == 0;
              boundary[j_dim][1] = row == n_div;
              block += row*2*strides_block[j_dim];
            }
            // there are `n_corner` possible surface elements that share vertex `block` with normal in direction `i_dim`
            int n_vert = math::pow(3, n_dim - 1); // these elements collectively have `n_vert` vertices
            std::vector<int> verts(n_vert); // `i_block` for each of said vertices
            // populate `verts` and compute search direction
            for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
              // compute block index
              int vert = block;
              for (int j_dim = 0, k_dim = 0; j_dim < n_dim; ++j_dim) if (j_dim != i_dim) {
                int row_incr = ((i_vert/math::pow(3, n_dim - 2 - k_dim))%3 - 1);
                vert += (!boundary[j_dim][row_incr > 0])*((i_vert/math::pow(3, n_dim - 2 - k_dim))%3 - 1)*strides_block[j_dim];
                ++k_dim;
              }
              // find index of vertex in `i_contour`
              int i_con;
              if (i_contour(vert) < 0) { // vertex isn't in contour yet, so add it
                i_con = i_block.size();
                i_block.push_back(vert);
                directions.push_back(Eigen::VectorXd::Zero(n_dim));
                i_contour(vert) = i_con;
              } else i_con = i_contour(vert); // vertex is already in contour, so just fetch its index
              verts[i_vert] = i_con;
              // update search direction so that it is not parallel to face
              directions[i_con][i_dim] += (1 - 2*(sample[sample0] > value));;
            }
            // does the orientation of this face need to be flipped to be consistent with the surface normal vector?
            bool flip = (sample[sample0] > value) != (i_dim == 1);
            // add faces
            for (int i_elem = 0; i_elem < n_corner; ++i_elem) {
              // add vertex indices
              for (int i_corner = 0; i_corner < n_corner; ++i_corner) {
                int i_vert = 0;
                for (int j_dim = 0; j_dim < n_dim - 1; ++j_dim) {
                  int stride = math::pow(2, n_dim - 2 - j_dim);
                  i_vert += ((i_elem/stride)%2 + (i_corner/stride)%2)*math::pow(3, n_dim - 2 - j_dim);
                }
                faces.push_back(verts[i_vert]);
              }
              // flip orientation if necessary
              if (flip) {
                auto start = faces.end() - n_corner;
                for (int i_row = 0; i_row < n_corner/2; ++i_row) {
                  std::swap(start[2*i_row], start[2*i_row + 1]);
                }
              }
            }
          }
        }
      }
    }
  }
  // compute initial reference coordinates of contour vertices, straight from the candidate block
  // these will need to be adjusted to lie exactly on the contour surface
  con.vert_ref_coords.resize(i_block.size(), n_dim);
  for (unsigned i = 0; i < i_block.size(); ++i) {
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      con.vert_ref_coords(i, i_dim) = ((i_block[i]/strides_block[i_dim])%(2*n_div + 1))/(2.*n_div);
    }
  }
  // compute gradient at quadrature points (used for projection)
  Eigen::VectorXd gradient(n_qpoint*n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    gradient(Eigen::seqN(i_dim*n_qpoint, n_qpoint)) = math::dimension_matvec(bas.diff_mat(), vars(Eigen::seqN((n_var - 1)*n_qpoint, n_qpoint)), i_dim);
  }
  // project points to contour surface
  // move in line search direction (computed above) and compute distance to move with newton's method
  for (int i_newton = 0; i_newton < n_newton; ++i_newton) {
    for (int i_vert = 0; i_vert < con.vert_ref_coords.rows(); ++i_vert) {
      auto coords = con.vert_ref_coords(i_vert, Eigen::all);
      double curr_value = sample_qpoint_data(vars(Eigen::seqN((n_var - 1)*n_qpoint, n_qpoint)), coords)(0);
      Eigen::VectorXd grad = sample_qpoint_data(gradient, coords).transpose(); // interpolate gradient to current coordinates
      auto dir = directions[i_vert].transpose();
      double diff = (value - curr_value)/(dir*grad + 1e-4*grad.norm());
      diff = std::max(-.5/n_div, std::min(.5/n_div, diff)); // limit search distance to prevent crazy-looking contours
      coords += dir*diff; // small stabilization term accounts for cases where both numerator and denominator -> 0
    }
  }
  // fetch jacobian (used for normal calculation)
  Eigen::VectorXd qpoint_jac(n_dim*n_dim*n_qpoint);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        qpoint_jac((i_dim*n_dim + j_dim)*n_qpoint + i_qpoint) = el.jacobian(i_dim, j_dim, i_qpoint);
      }
    }
  }
  // compute normals
  con.normals.resize(i_block.size(), n_dim);
  for (unsigned i_vert = 0; i_vert < i_block.size(); ++i_vert) {
    Eigen::MatrixXd coords = con.vert_ref_coords(i_vert, Eigen::all);
    Eigen::VectorXd grad = sample_qpoint_data(gradient, coords).transpose();
    Eigen::MatrixXd jac_t = sample_qpoint_data(qpoint_jac, coords); // 1 by n_dim*n_dim
    jac_t.resize(n_dim, n_dim); // automatically transposed bc of storage order
    con.normals(i_vert, Eigen::all) = (jac_t.householderQr().solve(grad)).normalized();
  }
  // put face info into Eigen matrix
  con.elem_vert_inds.resize(faces.size()/n_corner, n_corner);
  for (int i_elem = 0; i_elem < int(faces.size())/n_corner; ++i_elem) {
    for (int i_corner = 0; i_corner < n_corner; ++i_corner) {
      con.elem_vert_inds(i_elem, i_corner) = faces[i_elem*n_corner + i_corner];
    }
  }
  return con;
}
#endif

}
