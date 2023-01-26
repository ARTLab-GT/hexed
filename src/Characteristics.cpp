#include <Characteristics.hpp>

namespace hexed
{

double Characteristics::nrml(Mat<> vec)
{
  return dir.dot(vec);
}

Mat<> Characteristics::tang(Mat<> vec)
{
  return vec - dir*nrml(vec);
}

Characteristics::Characteristics(Mat<> state, Mat<> direction)
: dir{direction/direction.norm()},
  n_dim{int(direction.size())},
  mass{state(n_dim)},
  veloc{state(Eigen::seqN(0, n_dim))/mass}
{
  int n_dim = dir.size();
  // compute more properties of the reference state
  double vsq = veloc.squaredNorm();
  double pres = .4*(state(n_dim + 1) - .5*mass*vsq);
  double sound_speed = std::sqrt(1.4*pres/mass);
  // compute eigenvalues
  vals(2) = nrml(veloc);
  vals(0) = vals(2) - sound_speed;
  vals(1) = vals(2) + sound_speed;
  // compute 1D eigenvectors
  double d_mass = 1;
  for (int sign = 0; sign < 2; ++sign) {
    double d_veloc = (2*sign - 1)*sound_speed/mass*d_mass;
    double d_pres = 1.4*pres/mass*d_mass;
    vecs(Eigen::all, sign) <<
      d_mass*vals(2) + mass*d_veloc,
      d_mass,
      d_pres/.4 + .5*d_mass*vsq + mass*vals(2)*d_veloc;
  }
  vecs(Eigen::all, 2) << d_mass*vals(2), d_mass, .5*d_mass*vsq;
}

Mat<dyn, 3> Characteristics::decomp(Mat<> state)
{
  Mat<> mmtm = state(Eigen::seqN(0, n_dim));
  // component of tangential momentum perturbation which is not induced by mass perturbation
  Mat<> mmtm_correction = tang(mmtm) - state(n_dim)*tang(veloc);
  // compute state for 1D eigenvector problem
  Mat<3> state_1d;
  state_1d <<
    nrml(mmtm),
    state(n_dim),
    state(n_dim + 1) - veloc.dot(mmtm_correction);
  // decompose 1D state into eigenvectors
  Mat<1, 3> eig_basis = vecs.colPivHouseholderQr().solve(state_1d).transpose();
  Mat<3, 3> eig_decomp = vecs.array().rowwise()*eig_basis.array();
  // ND eigenvector decomposition
  Mat<dyn, 3> d(state.rows(), 3);
  d(Eigen::seqN(n_dim, 2), Eigen::all) = eig_decomp(Eigen::seqN(1, 2), Eigen::all);
  d(Eigen::seqN(0, n_dim), Eigen::all) = dir*eig_decomp(0, Eigen::all) + tang(veloc)*eig_basis; // second term accounts for tangential momentum induced by mass perturbation
  d(Eigen::seqN(0, n_dim), 2) += mmtm_correction;
  d(n_dim + 1, 2) += veloc.dot(mmtm_correction); // correct energy to account for tangential momentum perturbation
  return d;
}

}
