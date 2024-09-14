#include <hexed/Namespace.hpp>

namespace hexed {

bool Namespace::exists(std::string name) {
  return _ints.count(name) || _doubles.count(name) || _strings.count(name) || _arrays.count(name);
}

bool Namespace::exists_recursive(std::string name) {
  if (exists(name)) return true;
  if (!supers.empty()) {
    auto predicate = [name](std::shared_ptr<Namespace>& space) {return space->exists_recursive(name);};
    return std::all_of(supers.begin(), supers.end(), predicate);
  }
  return false;
}

std::vector<std::string> Namespace::names() const {
  std::vector<std::string> n;
  for (auto& pair : _ints)    n.push_back(pair.first);
  for (auto& pair : _doubles) n.push_back(pair.first);
  for (auto& pair : _strings) n.push_back(pair.first);
  for (auto& pair : _arrays) n.push_back(pair.first);
  std::sort(n.begin(), n.end());
  return n;
}

void Namespace::assign_array(Array<double> assign_to, std::string name) {
  auto arr = lookup<Array<double>>(name);
  if (arr) assign_to = *arr;
  else {
    auto d = lookup<double>(name);
    if (d) assign_to = *d;
    else HEXED_THROW(format_str(1000, "No numeric variable named `%s`.", name.c_str()));
  }
}

}
