#include <hexed/Deformed_element.hpp>
#include <hexed/math.hpp>
#include <hexed/Tree.hpp>

namespace hexed {

Deformed_element::Deformed_element(Storage_params params, Tree& t, int aniso_r_level)
: Element(params, t, true, aniso_r_level, true)
, n_qpoint{params.n_qpoint()}
, jac_dat{(n_dim*n_dim + 1)*n_qpoint}
{}

void Deformed_element::set_jacobian(const Basis& basis) {
  auto diff_mat = basis.diff_mat();
  const int n_qpoint = params.n_qpoint();
  // compute jacobian
  Eigen::VectorXd jac(n_dim*n_dim*n_qpoint);
  Array<double> shape_pos = position(basis);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      auto jac_entry {jac.segment((i_dim*n_dim + j_dim)*n_qpoint, n_qpoint)};
      jac_entry = math::dimension_matvec(diff_mat, shape_pos(i_dim).vector(), j_dim)/tree->nominal_shape()(j_dim);
    }
  }
  // compute interior normals
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    Eigen::MatrixXd qpoint_jac(n_dim, n_dim);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        qpoint_jac(i_dim, j_dim) = jac((i_dim*n_dim + j_dim)*n_qpoint + i_qpoint);
      }
    }
    jac_dat(n_dim*n_dim*n_qpoint + i_qpoint) = qpoint_jac.determinant();
    //if (!(jac_dat(n_dim*n_dim*n_qpoint + i_qpoint) > 0.)) printf("WARNING: nonpositive jacobian!\n");
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      Eigen::MatrixXd copy = qpoint_jac;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        copy(Eigen::all, i_dim).setUnit(j_dim);
        jac_dat((i_dim*n_dim + j_dim)*n_qpoint + i_qpoint) = copy.determinant();
      }
    }
  }
  // write surface normals to faces
  int nfq = n_qpoint/params.row_size;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      Eigen::MatrixXd bound_mat = basis.boundary()(sign, Eigen::all);
      Eigen::MatrixXd face_jac(nfq, n_dim*n_dim);
      for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac) {
        face_jac(Eigen::all, i_jac)
          = math::dimension_matvec(bound_mat, jac(Eigen::seqN(i_jac*n_qpoint, n_qpoint)), i_dim);
      }
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
        Eigen::MatrixXd qpoint_jac(n_dim, n_dim);
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          for (int k_dim = 0; k_dim < n_dim; ++k_dim) {
            qpoint_jac(j_dim, k_dim) = face_jac(i_qpoint, j_dim*n_dim + k_dim);
          }
        }
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          qpoint_jac(Eigen::all, i_dim).setUnit(j_dim);
          HEXED_ASSERT(face(2*i_dim + sign).normal().shape()[0] == n_dim, "normal has wrong size")
          face(2*i_dim + sign).normal()(j_dim)[i_qpoint] = qpoint_jac.determinant();
        }
      }
    }
  }
  // set local TSS
  auto bound_mat = basis.boundary();
  Eigen::MatrixXd vertex_nrml (params.n_vertices(), n_dim*n_dim);
  for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac) {
    // extrapolate one entry of the Jacobian to the vertex
    Eigen::Map<Eigen::VectorXd> rln(reference_level_normals() + i_jac*params.n_qpoint(), n_qpoint);
    vertex_nrml.col(i_jac) = math::hypercube_matvec(bound_mat, rln);
  }
  Eigen::Map<Eigen::VectorXd> jac_det(jacobian_determinant(), n_qpoint);
  Eigen::VectorXd vertex_det = math::hypercube_matvec(bound_mat, jac_det);
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    double norm_sum = 0.;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      double norm_sq = 0.;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        double coef = vertex_nrml(i_vert, i_dim*n_dim + j_dim);
        norm_sq += coef*coef;
      }
      norm_sum += std::sqrt(norm_sq)/nominal_shape(i_dim);
    }
    // for deformed elements this is a essentially a measure of the amount of stretching in each dimension
    vertex_time_step_scale(i_vert) = vertex_det(i_vert)/norm_sum;
  }
}

double* Deformed_element::reference_level_normals() {
  return jac_dat.data();
}

double* Deformed_element::jacobian_determinant() {
  return jac_dat.data() + n_dim*n_dim*n_qpoint;
}

double Deformed_element::jacobian(int i_dim, int j_dim, int i_qpoint) {
  Eigen::MatrixXd inv(n_dim, n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      inv(i_dim, j_dim) = jac_dat((i_dim*n_dim + j_dim)*n_qpoint + i_qpoint);
    }
  }
  return inv.inverse()(i_dim, j_dim)*jac_dat(n_dim*n_dim*n_qpoint + i_qpoint);
}

double Deformed_element::jacobian_determinant(int i_qpoint) {
  return jac_dat(n_dim*n_dim*n_qpoint + i_qpoint);
}

bool Deformed_element::deformed() const {return true;}

}
