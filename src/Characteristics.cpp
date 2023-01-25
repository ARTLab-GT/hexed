#include <Characteristics.hpp>
#include <iostream>

namespace hexed
{

Characteristics::Characteristics(Mat<> state, Mat<> direction)
: dir{direction/direction.norm()}, n_dim{int(direction.size())}
{
  int n_dim = dir.size();
  double mass = state(n_dim);
  auto mmtm = state(Eigen::seqN(0, n_dim));
  double pres = .4*(state(n_dim + 1) - .5*mmtm.squaredNorm()/mass);
  double sound_speed = std::sqrt(1.4*pres/mass);
  vals(2) = dir.dot(mmtm/mass);
  vals(0) = vals(2) - sound_speed;
  vals(1) = vals(2) + sound_speed;
  double d_mass = 1;
  for (int sign = 0; sign < 2; ++sign) {
    double d_veloc = (2*sign - 1)*sound_speed/mass*d_mass;
    double d_pres = 1.4*pres/mass*d_mass;
    vecs(Eigen::all, sign) <<
      d_mass*vals(2) + mass*d_veloc,
      d_mass,
      d_pres/.4 + .5*d_mass*vals(2)*vals(2) + mass*vals(2)*d_veloc;
  }
  vecs(Eigen::all, 2) << d_mass*vals(2), d_mass, .5*d_mass*vals(2)*vals(2);
  std::cout << vals << "\n\n" << vecs << "\n\n";
}

Mat<dyn, 3> Characteristics::decomp(Mat<> state)
{
  Mat<3> state_1d;
  state_1d <<
    dir.dot(state(Eigen::seqN(0, n_dim))),
    state(n_dim),
    state(n_dim + 1);
  Mat<1, 3> eig_basis = vecs.colPivHouseholderQr().solve(state_1d).transpose();
  Mat<3, 3> eig_decomp = vecs.array().rowwise()*eig_basis.array();
  Mat<dyn, 3> d(state.rows(), 3);
  d(Eigen::seqN(n_dim, 2), Eigen::all) = eig_decomp(Eigen::seqN(1, 2), Eigen::all);
  d(Eigen::seqN(0, n_dim), Eigen::all) = dir*eig_decomp(0, Eigen::all);
  d(Eigen::seqN(0, n_dim), 0) += state(Eigen::seqN(0, n_dim)) - dir*state_1d(0);
  std::cout << eig_basis << "\n\n" << eig_decomp << "\n\n" << d << "\n\n";
  return d;
}

}
