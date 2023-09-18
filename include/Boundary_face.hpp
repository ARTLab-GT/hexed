#ifndef HEXED_BOUNDARY_FACE_HPP_
#define HEXED_BOUNDARY_FACE_HPP_

#include <memory>
#include "math.hpp"
#include "Storage_params.hpp"
#include "Basis.hpp"

namespace hexed
{

/*!
 * Represents an element face where a boundary condition is to be applied without details
 * about the element or connection direction
 */
class Boundary_face
{
  public:
  virtual ~Boundary_face() = default;
  virtual Storage_params storage_params() = 0;
  virtual double* ghost_face() = 0;
  virtual double* inside_face() = 0;
  virtual int i_dim() = 0;
  virtual bool inside_face_sign() = 0;
  virtual double* surface_normal() = 0; //!< \note has to have a name that's different from `Face_connection`
  virtual double* surface_position() = 0;
  virtual double* state_cache() = 0; //!< \brief can be used to record the state for state-dependent flux boundary conditions
  virtual double* flux_cache() = 0; //!< \brief is used to record the flux for surface visualization/integrals
};

}
#endif
