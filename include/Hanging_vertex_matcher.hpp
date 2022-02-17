#ifndef HANGING_VERTEX_MATCHER_HPP_
#define HANGING_VERTEX_MATCHER_HPP_

#include "Element.hpp"

namespace cartdg
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

  public:
  // Note: sensitive to the order in which elements are provided
  Hanging_vertex_matcher(std::vector<Element*> fine_elements, int i_dim, bool is_positive);
  /*
   * Sets the values at the hanging vertices to the proper interpolated values. The interpolation
   * is based on the values at the non-hanging vertices of the fine elements, not the
   * `Vertex::shared_max_value`. So, if the shared value is desired, then this function should
   * be called after `Element::fetch_shareable_value`.
   */
  void match(Element::shareable_value_access);
};

}
#endif
