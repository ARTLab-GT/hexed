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
  r += format(double_format, mmtm_res);
  r += format(double_format, mass_res);
  r += format(double_format, ener_res);
  r += format(double_format, adv_res);
  r += format(double_format, diff_res);
  r += format(double_format, flow_time);
  r += format(double_format, time_step);
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

void Iteration_status::set_time()
{
  start = std::chrono::system_clock::now();
}

double Iteration_status::start_time()
{
  return std::chrono::duration_cast<std::chrono::nanoseconds>(start.time_since_epoch()).count()*1e-9;
}

double Iteration_status::wall_time()
{
  return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start).count()*1e-9;
}

}
