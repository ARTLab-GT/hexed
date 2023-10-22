#ifndef HEXED_ITERATIVE_HPP_
#define HEXED_ITERATIVE_HPP_

#include "Linear_equation.hpp"

//! \brief iterative solvers for linear equations, in particular Krylov subspace methods
namespace hexed::iterative
{

/*! \brief solves a linear equation approximately with
 * [GMRES](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method)
 * \details GMRES is restarted every `n_restart` iterations.
 * A total of `n_iters` restarts are performed, requiring approximately `n_restart*n_iters` matrix-vector products.
 * Requires the `Linear_equation` object to have storage for at least `n_restart + 3` vectors.
 */
void gmres(Linear_equation& equation, int n_restart, int n_iters);

/*! \brief solves a linear equation approximately with
 * [BiCGStab](https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method)
 * \details Runs `n_iters` iterations.
 * Requires the `Linear_equation` object to storage for at least 8 vectors.
 */
void bicgstab(Linear_equation& equation, int n_iters);

}
#endif
