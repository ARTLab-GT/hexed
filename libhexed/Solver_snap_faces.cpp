#include <Solver.hpp>
#include <Spatial.hpp>
#include <pde.hpp>
#include <kernels.hpp>
#include <Gauss_lobatto.hpp>

// this gets its own file cause it's its own special brand of scuffed
// for example, it invokes the Neighbor kernel direcly... honestly i should really just rewrite this function
void hexed::Solver::snap_faces()
{
  // perform basic snapping
  const int nfq = params.n_qpoint()/params.row_size;
  auto& bc_cons {acc_mesh->boundary_connections()};
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    Lock::Acquire acq(bc_cons[i_con].element().lock);
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh->boundary_condition(bc_sn).mesh_bc->snap_node_adj(bc_cons[i_con], basis);
  }
  // ...and that almost works by itself, but in 3D the presence of hanging nodes and ill-behaved warping correction
  // can result in faces matching up imperfectly.
  // So now, we have to coerce neighboring faces to match exactly, which is more of a pain than you would think
  if (params.n_dim < 3) return;
  // create a Gauss-lobatto basis and interpolation matrices between bases
  Gauss_lobatto lob(std::max(2, basis.row_size - 1));
  Mat<dyn, dyn> to_lob = basis.interpolate(lob.nodes());
  Mat<dyn, dyn> from_lob = lob.interpolate(basis.nodes());
  // basic setup
  _put_cache();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].state();
    for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) state[i_dof] = 0;
  }
  // Lift the node adjustments into the interior state so that we can use the kernels to operate on them.
  // The entire row of quadrature points that align with a given boundary quadrature point will have the same value of the node adjustment
  // written to the first state variable.
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    auto& con = bc_cons[i_con];
    Lock::Acquire acq(con.element().lock);
    int bc_sn = con.bound_cond_serial_n();
    if (bc_sn == acc_mesh->surface_bc_sn()) {
      double* state = con.element().state();
      double* node_adj = con.element().node_adjustments();
      if (node_adj) {
        node_adj += (2*con.i_dim() + con.inside_face_sign())*params.n_qpoint()/params.row_size;
        for (Row_index ind(params.n_dim, params.row_size, con.i_dim()); ind; ++ind) {
          for (int i_node = 0; i_node < params.row_size; ++i_node) {
            // In order to compare node adjustments for neighboring elements, the signs have to be consistent.
            // The `math::sign` term flips some of them so that a positive node adjustment always points out of the domain.
            state[ind.i_qpoint(i_node)] = math::sign(con.inside_face_sign())*node_adj[ind.i_face_qpoint()];
          }
        }
      }
    }
  }
  // write to the faces
  compute_write_face(_kernel_mesh());
  compute_prolong(_kernel_mesh());
  // In the case of the hanging node connections, we want both sides to use the coarse face's value,
  // since the coarse element cannot exactly represent values from the fine element.
  // Here we mess with the values on both sides so that the average is equal to the coarse value.
  auto& ref_cons = acc_mesh->deformed().refined_connections();
  #pragma omp parallel for
  for (int i_con = 0; i_con < ref_cons.size(); ++i_con) {
    for (int i_fine = 0; i_fine < 2; ++i_fine) {
      double* state = ref_cons[i_con].connection(i_fine).state(0, false);
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) state[i_qpoint] *= 2;
      state = ref_cons[i_con].connection(i_fine).state(1, false);
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) state[i_qpoint] *= 0;
    }
  }
  // just to be safe, make sure the boundary faces have the right values
  Copy fake_bc;
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    auto& con = bc_cons[i_con];
    Lock::Acquire acq(con.element().lock);
    fake_bc.apply_state(con);
  }
  // trick the neighbor kernel into taking the averages for us
  (*kernel_factory<Spatial<pde::Smooth_art_visc, true>::Neighbor>(params.n_dim, params.row_size, 0, 0, 1., 1.))(acc_mesh->deformed().kernel_connections());
  compute_restrict(_kernel_mesh(), false, true);
  // At this point, the LDG face state has the averaged values.
  // Write the original values to the regular face state for comparison.
  compute_write_face(_kernel_mesh());
  // finally, make the corrections to the node adjustments
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    auto& con = bc_cons[i_con];
    Lock::Acquire acq(con.element().lock);
    int bc_sn = con.bound_cond_serial_n();
    if (bc_sn == acc_mesh->surface_bc_sn()) {
      double* state = con.element().state();
      double* node_adj = con.element().node_adjustments();
      if (node_adj) {
        node_adj += (2*con.i_dim() + con.inside_face_sign())*params.n_qpoint()/params.row_size;
        for (int j_dim = 0; j_dim < params.n_dim; ++j_dim) if (j_dim != con.i_dim()) {
          for (Row_index ind(params.n_dim, params.row_size, j_dim); ind; ++ind) {
            // obtain the corrections from the faces and apply them to the node adjustments stored in the interior state
            Mat<> row(params.row_size);
            for (int i_node = 0; i_node < params.row_size; ++i_node) {
              row(i_node) = state[ind.i_qpoint(i_node)];
            }
            Mat<> lrow = to_lob*row;
            for (int i_sign = 0; i_sign < 2; ++i_sign) {
              // use the difference demanded by each face instead of the value demanded by each face,
              // since that avoids a problem with hanging nodes
              double diff =   con.element().face(2*j_dim + i_sign, true)[ind.i_face_qpoint()]
                            - con.element().face(2*j_dim + i_sign, false)[ind.i_face_qpoint()];
              lrow(i_sign*(lob.row_size - 1)) += diff;
            }
            row = from_lob*lrow;
            for (int i_node = 0; i_node < params.row_size; ++i_node) {
              state[ind.i_qpoint(i_node)] = row(i_node);
            }
          }
        }
        // transfer values from interior state to the boundary faces where it belongs
        for (Row_index ind(params.n_dim, params.row_size, con.i_dim()); ind; ++ind) {
          node_adj[ind.i_face_qpoint()] = math::sign(con.inside_face_sign())*state[ind.i_qpoint(0)];
        }
      }
    }
  }
  // clean up
  _get_cache();
}

