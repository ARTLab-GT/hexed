#ifndef HEXED_ITERATIVE_HPP_
#define HEXED_ITERATIVE_HPP_

#include "Linear_equation.hpp"

namespace hexed::iterative
{

void gmres(Linear_equation& equation, int n_restart, int n_iters);
void bicgstab(Linear_equation& equation, int n_iters);

}
#endif
