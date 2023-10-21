#include <Solver.hpp>
#include <Prolong_refined.hpp>

namespace hexed
{

Solver::Linearized::Linearized(Solver& s)
: _solver{s}
{
  Mat<> _weights = math::pow_outer(_solver.basis.node_weights(), _solver.params.n_dim);
  Mat<> _ref_state(_solver.params.n_var);
  double mmtm = 0.;
  for (int i_dim = 0; i_dim < _solver.params.n_dim; ++i_dim) {
    double component = _solver._namespace->lookup<double>("freestream" + std::to_string(_solver.params.n_dim)).value();
    mmtm += math::pow(component, 2);
  }
  mmtm = std::sqrt(mmtm);
  for (int i_dim = 0; i_dim < _solver.params.n_dim; ++i_dim) _ref_state(i_dim) = mmtm;
  for (int i_var = _solver.params.n_dim; i_var < _solver.params.n_var; ++i_var) {
    _ref_state(i_var) = _solver._namespace->lookup<double>("freestream" + std::to_string(i_var)).value();
  }
}

int Solver::Linearized::n_vecs() {return n_storage;}

void Solver::Linearized::scale(int output, int input, double scalar)
{
  auto& elems = _solver.acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* out = elems[i_elem].stage(storage_start + output);
    double* in  = elems[i_elem].stage(storage_start + input);
    for (int i_dof = 0; i_dof < _solver.params.n_dof(); ++i_dof) out[i_dof] = scalar*in[i_dof];
  }
}

void Solver::Linearized::add(int output, double coef0, int vec0, double coef1, int vec1)
{
  auto& elems = _solver.acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* out = elems[i_elem].stage(storage_start + output);
    double* in0 = elems[i_elem].stage(storage_start + vec0);
    double* in1 = elems[i_elem].stage(storage_start + vec1);
    for (int i_dof = 0; i_dof < _solver.params.n_dof(); ++i_dof) {
      out[i_dof] = coef0*in0[i_dof] + coef1*in1[i_dof];
    }
  }
}

double Solver::Linearized::inner(int input0, int input1)
{
  auto& elems = _solver.acc_mesh.elements();
  double total = 0;
  const int nq = _solver.params.n_qpoint();
  #pragma omp parallel for reduction (+:total)
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* in0 = elems[i_elem].stage(storage_start + input0);
    double* in1 = elems[i_elem].stage(storage_start + input1);
    for (int i_var = 0; i_var < _solver.params.n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        int i_dof = i_var*nq + i_qpoint;
        total += in0[i_dof]*_weights(i_qpoint)/_ref_state(i_var)*in1[i_dof];
      }
    }
  }
  return total;
}

void Solver::Linearized::matvec(int output, int input)
{
  std::unique_ptr<Kernel<Refined_face&>> prolong(kernel_factory<Prolong_refined>(
    _solver.params.n_dim, _solver.params.row_size, _solver.basis));
  add(-storage_start, 1., -storage_start, finite_diff, input);
  (*_solver.write_face)(_solver.acc_mesh.elements());
  (*prolong)(_solver.acc_mesh.refined_faces());
  _solver.update();
  add(output, 1/finite_diff, -storage_start, -1/finite_diff, -storage_start + 3);
  scale(-storage_start, -storage_start + 3, 1.);
  (*_solver.write_face)(_solver.acc_mesh.elements());
  (*prolong)(_solver.acc_mesh.refined_faces());
}

}
