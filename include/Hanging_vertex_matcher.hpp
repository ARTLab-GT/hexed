#ifndef HEXED_HANGING_VERTEX_MATCHER_HPP_
#define HEXED_HANGING_VERTEX_MATCHER_HPP_

#include "Element.hpp"

namespace hexed
{

/*
 * When hanging nodes are present, data at the hanging vertices (a.k.a. "hanging nodes")
 * must match data on the coarse face, interpolated to the location of the hanging vertices.
 * This class ensures that is the case by performing the necessary interpolation.
 */
class Hanging_vertex_matcher
{
  std::vector<Element*> elements;
  int id;
  bool isp;
  int n_dim;
  int n_vert;
  std::array<bool, 2> str;

  public:
  // Note: sensitive to the order in which elements are provided
  Hanging_vertex_matcher(std::vector<Element*> fine_elements, int i_dim, bool is_positive, std::array<bool, 2> stretch = {false, false});
  /*
   * Sets the values at the hanging vertices to the proper interpolated values.
   * Note: the interpolation is based on the values at the non-hanging vertices
   * of the fine elements, not the `Vertex::shared_max_value`.
   * So, if the shared value is desired, then this function should
   * be called after `Element::fetch_shareable_value`.
   */
  void match(Element::vertex_value_access);
};

}
#endif
