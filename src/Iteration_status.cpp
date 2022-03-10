#include <Iteration_status.hpp>

namespace cartdg
{

int Iteration_status::width()
{
  int w = number_width;
  for (auto label : labels) w = std::max<int>(w, label.size());
  return w;
}

std::string Iteration_status::header()
{
  std::string h = "";
  for (auto label : labels) {
    h += format("s", label.c_str());
  }
  for (int i = 0; i < n_sep; ++i) h.pop_back(); // delete trailing spaces
  return h;
}

std::string Iteration_status::report()
{
  std::string r = "";
  r += format("i", iteration);
  r += format(".15f", time);
  r += format(".15f", time_step);
  for (int i = 0; i < n_sep; ++i) r.pop_back();
  return r;
}

}
