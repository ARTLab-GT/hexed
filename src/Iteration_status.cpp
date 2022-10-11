#include <Iteration_status.hpp>

namespace hexed
{

int Iteration_status::width()
{
  int w = number_width;
  for (auto label : labels) w = std::max<int>(w, label.size());
  return w;
}

std::string Iteration_status::value_string()
{
  std::string r = "";
  r += format("i", iteration);
  r += format(".14e", flow_time);
  r += format(".14e", time_step);
  r += format("i", fix_admis_iters);
  return r;
}

std::string Iteration_status::header()
{
  std::string h = "";
  for (auto label : labels) {
    h += format("s", label.c_str());
  }
  for (unsigned i = 0; i < sep.size(); ++i) h.pop_back(); // delete trailing separator
  return h;
}

std::string Iteration_status::report()
{
  auto r = value_string();
  for (unsigned i = 0; i < sep.size(); ++i) r.pop_back();
  return r;
}

}
