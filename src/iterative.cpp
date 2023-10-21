#include <iterative.hpp>

namespace hexed::iterative
{

void gmres(Linear_equation& equation, int n_restart, int n_iters)
{
  const int krylov_start = 3;
  HEXED_ASSERT(equation.n_vecs() >= krylov_start + n_restart, "`Linear_equation` object has insufficient storage for the GMRES algorithm");
  for (int restart = 0; restart < n_iters; ++restart) {
    equation.matvec(2, 0);
    equation.add(2, -1., 2, 1., 1);
    double res_norm = equation.norm(2);
    equation.scale(krylov_start, 2, 1./res_norm);
    Mat<dyn, dyn> hessenberg = Mat<dyn, dyn>::Zero(n_restart + 1, n_restart);
    for (int iter = 0; iter < n_restart; ++iter) {
      int next = krylov_start + iter + 1;
      equation.matvec(next, krylov_start + iter);
      for (int col = 0; col < iter + 1; ++col) {
        hessenberg(col, iter) = equation.inner(next, krylov_start + col);
        equation.add(next, 1., next, -hessenberg(col, iter), krylov_start + col);
      }
      hessenberg(iter + 1, iter) = equation.norm(next);
      equation.scale(next, next, 1./hessenberg(iter + 1, iter));
    }
    Mat<> coefs = hessenberg.fullPivHouseholderQr().solve(res_norm*Mat<>::Unit(n_restart + 1, 0));
    for (int iter = 0; iter < n_restart; ++iter) {
      equation.add(0, 1., 0, coefs(iter), krylov_start + iter);
    }
  }
}

void bicgstab(Linear_equation& equation, int n_iters)
{
}

}
