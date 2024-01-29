#include <Solver.hpp>
#include <Spatial.hpp>
#include <pde.hpp>
#include <kernels.hpp>
#include <Gauss_lobatto.hpp>

// this gets its own file cause it's its own special brand of scuffed
// for example, it invokes the Neighbor kernel direcly... honestly i should really just rewrite this function
void hexed::Solver::snap_faces()
{
  const int nd = params.n_dim;
  const int nfq = params.n_qpoint()/params.row_size;
  auto& bc_cons {acc_mesh.boundary_connections()};
  auto& elems = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    Lock::Acquire acq(bc_cons[i_con].element().lock);
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).mesh_bc->snap_node_adj(bc_cons[i_con], basis);
  }
  Gauss_lobatto lob(std::max(2, basis.row_size - 1));
  Mat<dyn, dyn> to_lob = basis.interpolate(lob.nodes());
  Mat<dyn, dyn> from_lob = lob.interpolate(basis.nodes());
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].state();
    double* rk_ref = elems[i_elem].residual_cache();
    for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) rk_ref[i_dof] = state[i_dof];
  }
  for (int sweep = 0; sweep < nd - 1; ++sweep) {
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      double* state = elems[i_elem].state();
      for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) state[i_dof] = 0;
    }
    #pragma omp parallel for
    for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
      auto& con = bc_cons[i_con];
      Lock::Acquire acq(con.element().lock);
      int bc_sn = con.bound_cond_serial_n();
      if (bc_sn == acc_mesh.surface_bc_sn()) {
        double* state = con.element().state();
        double* node_adj = con.element().node_adjustments();
        if (node_adj) {
          node_adj += (2*con.i_dim() + con.inside_face_sign())*params.n_qpoint()/params.row_size;
          for (Row_index ind(params.n_dim, params.row_size, con.i_dim()); ind; ++ind) {
            for (int i_node = 0; i_node < params.row_size; ++i_node) {
              state[ind.i_qpoint(i_node)] = math::sign(con.inside_face_sign())*node_adj[ind.i_face_qpoint()];
            }
          }
        }
      }
    }
    compute_write_face(_kernel_mesh);
    compute_prolong(_kernel_mesh);
    Copy fake_bc;
    #pragma omp parallel for
    for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
      auto& con = bc_cons[i_con];
      Lock::Acquire acq(con.element().lock);
      fake_bc.apply_state(con);
    }
    if (sweep) {
      auto& ref_cons = acc_mesh.deformed().refined_connections();
      #pragma omp parallel for
      for (int i_con = 0; i_con < ref_cons.size(); ++i_con) {
        for (int i_fine = 0; i_fine < 2; ++i_fine) {
          double* state = ref_cons[i_con].connection(i_fine).state(0, false);
          for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) state[i_qpoint] *= 2;
          state = ref_cons[i_con].connection(i_fine).state(1, false);
          for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) state[nfq*params.n_var + i_qpoint] = 0;
        }
      }
    }
    (*kernel_factory<Spatial<pde::Smooth_art_visc, false>::Neighbor>(params.n_dim, params.row_size, 0, 1., 1.))(acc_mesh.deformed().kernel_connections());
    compute_restrict(_kernel_mesh, false, true);
    #pragma omp parallel for
    for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
      auto& con = bc_cons[i_con];
      Lock::Acquire acq(con.element().lock);
      int bc_sn = con.bound_cond_serial_n();
      if (bc_sn == acc_mesh.surface_bc_sn()) {
        double* state = con.element().state();
        double* node_adj = con.element().node_adjustments();
        if (node_adj) {
          node_adj += (2*con.i_dim() + con.inside_face_sign())*params.n_qpoint()/params.row_size;
          for (int j_dim = 0; j_dim < params.n_dim; ++j_dim) if (j_dim != con.i_dim()) {
            for (Row_index ind(params.n_dim, params.row_size, j_dim); ind; ++ind) {
              Mat<> row(params.row_size);
              for (int i_node = 0; i_node < params.row_size; ++i_node) {
                row(i_node) = state[ind.i_qpoint(i_node)];
              }
              Mat<> lrow = to_lob*row;
              for (int i_sign = 0; i_sign < 2; ++i_sign) {
                lrow(i_sign*(lob.row_size - 1)) = con.element().face(2*j_dim + i_sign, true)[ind.i_face_qpoint()];
              }
              row = from_lob*lrow;
              for (int i_node = 0; i_node < params.row_size; ++i_node) {
                state[ind.i_qpoint(i_node)] = row(i_node);
              }
            }
          }
          for (Row_index ind(params.n_dim, params.row_size, con.i_dim()); ind; ++ind) {
            node_adj[ind.i_face_qpoint()] = math::sign(con.inside_face_sign())*state[ind.i_qpoint(0)];
          }
        }
      }
    }
  }
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].state();
    double* rk_ref = elems[i_elem].residual_cache();
    for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) state[i_dof] = rk_ref[i_dof];
  }
}

