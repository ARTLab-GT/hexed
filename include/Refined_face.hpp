#ifndef CARTDG_REFINED_FACE_HPP_
#define CARTDG_REFINED_FACE_HPP_

#include <Storage_params.hpp>

namespace cartdg
{

class Refined_face
{
  protected:
  Eigen::VectorXd fine;
  double* coarse;

  public:
  Refined_face(Storage_params, double* coarse_face);
  double* fine_face(int i_face);
  double* coarse_face();
};

typedef std::vector<std::unique_ptr<Refined_face>> ref_face_vec;

}
#endif
