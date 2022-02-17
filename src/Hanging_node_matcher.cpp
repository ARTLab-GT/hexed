#include <Hanging_node_matcher.hpp>

namespace cartdg
{

Hanging_node_matcher::Hanging_node_matcher(std::vector<Element*> fine_elements, int i_dim, bool is_positive)
: elements{fine_elements}, id{i_dim}, isp{is_positive}
{}

void Hanging_node_matcher::match(Element::shareable_value_access access_func)
{
}

}
