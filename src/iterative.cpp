#include <iterative.hpp>

namespace hexed::iterative
{

void gmres(Linear_equation& equation, int n_restart, int n_iters)
{
  double res_norm = equation.norm(1);
  const int krylov_start = 3;
  equation.scale(krylov_start, 1, 1./res_norm);
  Mat<dyn, dyn> hessenberg = Mat<dyn, dyn>::Zero(n_restart + 1, n_restart);
  for (int iter = 0; iter < n_restart; ++iter) {
    equation.matvec(2, krylov_start + iter);
    for (int col = 0; col < iter + 1; ++col) {
      hessenberg(col, iter) = equation.inner(2, krylov_start + col);
      equation.add(2, 1., 2., -hessenberg(col, iter), krylov_start + col);
    }
    hessenberg(iter + 1, iter) = equation.norm(2);
    equation.scale(krylov_start + iter + 1, 2, 1./hessenberg(iter + 1, iter));
  }
  Mat<> coefs = hessenberg.fullPivHouseholderQr().solve(res_norm*Mat<>::Unit(n_restart + 1, 0));
  for (int iter = 0; iter < n_restart; ++iter) {
    equation.add(0, 1., 0, coefs(iter), krylov_start + iter);
  }
}

}
