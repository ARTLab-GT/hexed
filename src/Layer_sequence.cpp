#include <cmath>
#include <Layer_sequence.hpp>

namespace hexed
{

Layer_sequence::Layer_sequence(double wall_spacing, int n_wall)
: ws{wall_spacing},
  nw{n_wall},
  n_grow{int(std::ceil(std::log(.5/ws)/std::log(2)))}
{}

int Layer_sequence::n_layers()
{
  return nw + n_grow;
}

double Layer_sequence::spacing(int i_layer)
{
  double first_growth = std::pow(2., nw - n_layers());
  return (i_layer < nw) ? ws : std::pow(2., i_layer - n_layers())*(1 - ws*nw + first_growth);
}

double Layer_sequence::cumulative_height(int i_layer)
{
  // might as well just do it the brute force way
  double total = 0;
  for (int j_layer = 0; j_layer < i_layer; ++j_layer) total += spacing(j_layer);
  return total;
}

}
