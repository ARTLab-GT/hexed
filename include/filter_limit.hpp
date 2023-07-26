#ifndef FILTER_LIMIT_HPP_
#define FILTER_LIMIT_HPP_

#include "Basis.hpp"

namespace hexed
{

/*! \brief Applies a low-pass filter to polynomial data.
 * \details Limits the Legendre modes of multidimensional polynomial
 * such that an \f$n\f$th degree mode will have norm no more than `decay_rate`\f$^n\f$.
 * Each mode which violates this condition is scaled to the largest value that satisfies it.
 * Thus, if `decay_rate == 1` then the data will never be altered,
 * whereas if `decay_rate == 0` then all the modes except the 0th-degree will be annhilated.
 * `data` is assumed to point to the values of a scalar polynomial at the `n_dim`-dimensional
 * quadrature points associated with the provided basis.
 */
void filter_limit(int n_dim, int row_size, double* data, const Basis&, double decay_rate);

}
#endif
