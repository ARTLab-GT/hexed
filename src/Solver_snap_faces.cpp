#include <Solver.hpp>
#include <Gauss_lobatto.hpp>

// this gets its own file cause it's its own special brand of scuffed
// for example, it invokes the Neighbor kernel direcly... honestly i should really just rewrite this function
void hexed::Solver::snap_faces()
{
  const int nd = params.n_dim;
  const int nfq = params.n_qpoint()/params.row_size;
  auto& bc_cons {acc_mesh->boundary_connections()};
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    Lock::Acquire acq(bc_cons[i_con].element().lock);
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh->boundary_condition(bc_sn).mesh_bc->snap_node_adj(bc_cons[i_con], basis);
  }
  Gauss_lobatto lob(std::max(2, basis.row_size - 1));
  Mat<dyn, dyn> to_lob = basis.interpolate(lob.nodes());
  Mat<dyn, dyn> from_lob = lob.interpolate(basis.nodes());
  auto apply_mat = [&](Mat<dyn, dyn> mat) {
    auto ext_cons = acc_mesh->extruded_connections();
    #pragma omp parallel for
    for (int i_con = 0; i_con < ext_cons.size(); ++i_con) {
      auto& con = ext_cons[i_con];
      auto& elem = con.element(0);
      HEXED_ASSERT(!elem.tree, "oops this was supposed to be extruded");
      auto dir = con.get_direction();
      Eigen::Map<Mat<>> node_adj(elem.node_adjustments() + (2*dir.i_dim[0] + dir.face_sign[0])*nfq, nfq);
      node_adj = math::hypercube_matvec(mat, node_adj);
    }
  };
  apply_mat(basis.interpolate(lob.nodes()));
  apply_mat(lob.interpolate(basis.nodes()));
}

