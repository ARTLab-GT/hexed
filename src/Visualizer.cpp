#include <Visualizer.hpp>
#include <Xdmf_wrapper.hpp>
#include <Tecplot_file.hpp>
#include <Csv.hpp>

namespace hexed
{

std::unique_ptr<Visualizer> Visualizer::create(std::string format, int n_dim_geom, int n_dim_topo, std::string file_name,
                                               std::vector<std::string> variable_names, double time, elem_type elem_t)
{
  std::unique_ptr<Visualizer> visualizer;
  if (format == "xdmf") {
    #if HEXED_USE_XDMF
    visualizer.reset(new Xdmf_wrapper(n_dim_geom, n_dim_topo, file_name, variable_names, time, elem_t));
    #else
    HEXED_ASSERT(false, "`format = xdmf` requires `USE_XDMF ON`");
    #endif
  } else if (format == "tecplot") {
    #if HEXED_USE_TECPLOT
    visualizer.reset(new Tecplot_file(file_name, n_dim_geom, n_dim_topo, variable_names, time));
    #else
    HEXED_ASSERT(false, "`format = tecplot` requires `USE_TECPLOT ON`");
    #endif
  } else if (format == "csv") {
    std::vector<std::string> all_names;
    for (int i_dim = 0; i_dim < n_dim_geom; ++i_dim) all_names.push_back("pos" + std::to_string(i_dim));
    all_names.insert(all_names.end(), variable_names.begin(), variable_names.end());
    visualizer.reset(new Csv(file_name, all_names));
  } else HEXED_ASSERT(false, format_str(1000, "visualization format `%s` not recognized", format.c_str()));
  return visualizer;
}

std::unique_ptr<Visualizer> Visualizer::create(std::string format, int n_dim_geom, int n_dim_topo,
                                               std::string file_name, const Output_data& dat, double time, elem_type elem_t)
{
  int n_var = dat.n_var(n_dim_geom);
  std::vector<std::string> variable_names;
  for (int i_var = 0; i_var < n_var; ++i_var) variable_names.push_back(dat.variable_name(n_dim_geom, i_var));
  return create(format, n_dim_geom, n_dim_topo, file_name, variable_names, time, elem_t);
}

}
