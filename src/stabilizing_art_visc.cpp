#include <stabilizing_art_visc.hpp>
#include <Row_rw.hpp>
#include <constants.hpp>

namespace hexed
{

template <int n_dim, int row_size>
class Stab_art_visc : public Kernel<Kernel_element&>
{
  static constexpr int n_qpoint = math::pow(row_size, n_dim);
  double _char_speed;
  double _ramp_center;
  double _half_width = 0.5;
  Mat<row_size> _row_weights;
  Mat<n_qpoint> _qpoint_weights;
  Mat<n_qpoint/row_size> _face_weights;
  Mat<row_size> _proj;

  public:
  Stab_art_visc(const Basis& basis, double char_speed) :
    _char_speed{char_speed},
    _ramp_center{-4.25*std::log(row_size - 1)/std::log(10)},
    _row_weights{basis.node_weights()},
    _qpoint_weights{math::pow_outer(_row_weights, n_dim)},
    _face_weights{math::pow_outer(_row_weights, n_dim - 1)},
    _proj{basis.orthogonal(row_size - 1).cwiseProduct(_row_weights)}
  {}

  void operator()(Sequence<Kernel_element&>& elements)
  {
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
    {
      double* state = elements[i_elem].state();
      double indicator_vars [2][n_qpoint];
      double norm_sq [2] {};
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        indicator_vars[0][i_qpoint] = 1./state[n_dim*n_qpoint + i_qpoint];
        double mmtm_sq = 0;
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) mmtm_sq += math::pow(state[i_dim*n_qpoint + i_qpoint], 2);
        indicator_vars[1][i_qpoint] = 1./(state[(n_dim + 1)*n_qpoint + i_qpoint] - .5*mmtm_sq*indicator_vars[0][i_qpoint]);
        for (int i_var = 0; i_var < 2; ++i_var) {
          norm_sq[i_var] += math::pow(indicator_vars[i_var][i_qpoint], 2)*_qpoint_weights(i_qpoint);
        }
      }
      Mat<2> nonsmooth = Mat<2>::Zero();
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
          Mat<row_size, 2> row = Row_rw<2, row_size>::read_row(indicator_vars[0], ind);
          Mat<2> dot = row.transpose()*_proj;
          nonsmooth += dot.cwiseProduct(dot)*_face_weights(ind.i_face_qpoint());
        }
      }
      double indicator = 0;
      for (int i_var = 0; i_var < 2; ++i_var) indicator += nonsmooth(i_var)/(norm_sq[i_var]*n_dim);
      indicator = std::log(indicator)/std::log(10);
      if (indicator <= _ramp_center - _half_width) indicator = 0;
      else if (indicator < _ramp_center + _half_width) indicator = .5*(1 + std::sin(constants::pi*(indicator - _ramp_center)/2/_half_width));
      else indicator = 1;
      elements[i_elem].uncert() = (row_size - 1)*_char_speed*elements[i_elem].nominal_size()*indicator;
    }
  }
};

void stabilizing_art_visc(Kernel_mesh mesh, double char_speed)
{
  (*kernel_factory<Stab_art_visc>(mesh.n_dim, mesh.row_size, mesh.basis, char_speed))(mesh.elems);
}

}
