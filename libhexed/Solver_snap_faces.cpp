#include <Solver.hpp>
#include <Spatial.hpp>
#include <pde.hpp>
#include <kernels.hpp>
#include <Gauss_lobatto.hpp>

// this gets its own file cause it's its own special brand of scuffed
// for example, it invokes the Neighbor kernel direcly... honestly i should really just rewrite this function
void hexed::Solver::snap_faces()
{
}

