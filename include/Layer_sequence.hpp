#ifndef HEXED_LAYER_SEQUENCE_HPP_
#define HEXED_LAYER_SEQUENCE_HPP_

namespace hexed
{

class Layer_sequence
{
  public:
  Layer_sequence(double wall_spacing, int n_wall);
  int n_layers();
  double spacing(int i_layer);
  double cumulative_height(int i_layer);
};

}
#endif
