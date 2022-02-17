#include <Hanging_node_matcher.hpp>
#include <math.hpp>

namespace cartdg
{

Hanging_node_matcher::Hanging_node_matcher(std::vector<Element*> fine_elements, int i_dim, bool is_positive)
: elements{fine_elements}, id{i_dim}, isp{is_positive}, n_vert{elements[0]->storage_params().n_vertices()/2}
{}

void Hanging_node_matcher::match(Element::shareable_value_access access_func)
{
}

}
