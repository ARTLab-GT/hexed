#ifndef CARTDG_REFINED_FACE_HPP_
#define CARTDG_REFINED_FACE_HPP_

#include <memory>
#include <Eigen/Dense>
#include "Storage_params.hpp"

namespace cartdg
{

/*
 * Stores "mortar" data for connecting faces of different refinement levels (connections with hanging
 * nodes). Contains a pointer to the face data for the coarse element and storage to contain its
 * interpolation onto the fine element faces.
 */
class Refined_face
{
  protected:
  int face_size;
  int n_face;
  Eigen::VectorXd fine;
  double* coarse;

  public:
  Refined_face(Storage_params, double* coarse_face);
  virtual ~Refined_face() = default;
  // The following return pointers to data for a single face.
  double* fine_face(int i_face); // pointer to the mortar (owned by `this`) fine face storage
  double* coarse_face(); // pointer to the external (owned by element) coarse face storage
};

}
#endif
