#ifndef CARTDG_REFINED_FACE_HPP_
#define CARTDG_REFINED_FACE_HPP_

#include <memory>
#include <Eigen/Dense>
#include "Storage_params.hpp"

namespace cartdg
{

/*
 * Stores "mortar" data for connecting faces of different refinement levels (connections with hanging
 * nodes). To use, first construct a `Refined_face` object from the coarse face, then use `elem_con`
 * objects to connect the faces of the fine elements to the pointers returned by `fine_face`. Use the
 * prolong kernel to copy the polynomial in `coarse_face()` to the `fine_face`s, apply the neighbor
 * kernel as usual, and then use the restrict kernel to take the flux from the fine faces back to the
 * coarse face.
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

typedef std::vector<std::vector<std::unique_ptr<Refined_face>>> ref_face_vec;

}
#endif
