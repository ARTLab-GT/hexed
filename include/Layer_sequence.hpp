#ifndef HEXED_LAYER_SEQUENCE_HPP_
#define HEXED_LAYER_SEQUENCE_HPP_

namespace hexed
{

/*
 * Selects a growth profile for anisotropic wall layers.
 * That is, given a few high-level parameters,
 * it can compute the height that each layer should have for effective boundary layer resolution.
 */
class Layer_sequence
{
  double ws;
  int nw;
  int n_grow; // number of growth layers (after same-height wall layers)
  public:
  /*
   * Will generate `n_wall` layers with thickness `wall_spacing`
   * and then add exponentially growing layers on top until it reaches a total height of 1.
   * Uses the minimum number of layers required to keep the growth ratio about 2.
   */
  Layer_sequence(double wall_spacing, int n_wall);
  int n_layers();
  double spacing(int i_layer);
  double cumulative_height(int i_layer);
};

}
#endif
