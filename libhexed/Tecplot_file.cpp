#include <config.hpp>
#if HEXED_USE_TECPLOT
#include <TECIO.h>
#include <MASTER.h>
#include <Tecplot_file.hpp>
#include <math.hpp>

namespace hexed
{

Tecplot_file::Tecplot_file(std::string file_name, int n_dim, int n_dim_topo_arg, std::vector<std::string> variable_names, double time, double heat_rat, double gas_const)
: n_dim{n_dim}, n_dim_topo{n_dim_topo_arg}, n_var{int(variable_names.size())}, time{time}, strand_id{1}, i_zone{0}, file_handle{nullptr}
{
  std::string var_name_list = "";
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    var_name_list += " position" + std::to_string(i_dim);
  }
  for (int i_var = 0; i_var < n_var; ++i_var) {
    var_name_list += " " + variable_names[i_var];
  }
  var_name_list.erase(0, 1); // erase leading space
  tecFileWriterOpen(file_name.c_str(), "flow solution", var_name_list.c_str(), FILEFORMAT_SZL, FILETYPE_FULL, 0, nullptr, &file_handle);
  tecDataSetAddAuxData(file_handle, "Common.Incompressible", "False");
  tecDataSetAddAuxData(file_handle, "Common.Incompressible", "False");
  tecDataSetAddAuxData(file_handle, "Common.VectorVarsAreVelocity", "False");
  tecDataSetAddAuxData(file_handle, "Common.Gamma", std::to_string(heat_rat).c_str());
  tecDataSetAddAuxData(file_handle, "Common.GasConstant", std::to_string(gas_const).c_str());
  std::vector<std::string> xyz {"X", "Y", "Z"};
  std::vector<std::string> uvw {"U", "V", "W"};
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    tecDataSetAddAuxData(file_handle, ("Common." + xyz[i_dim] + "Var").c_str(), std::to_string(i_dim + 1).c_str());
    tecDataSetAddAuxData(file_handle, ("Common." + uvw[i_dim] + "Var").c_str(), std::to_string(n_dim + i_dim + 1).c_str());
  }
  tecDataSetAddAuxData(file_handle, "Common.CVar", std::to_string(2*n_dim + 1).c_str());
  tecDataSetAddAuxData(file_handle, "Common.DensityVar", std::to_string(2*n_dim + 1).c_str());
  tecDataSetAddAuxData(file_handle, "Common.StagnationEnergyVar", std::to_string(2*n_dim + 2).c_str());
}

void Tecplot_file::write_block(Array<double> pos, Array<double> vars)
{
  HEXED_ASSERT(pos(0).same_shape(vars(0)), "`pos` and `vars` must have the same shape");
  Structured_block block(*this, pos.shape()[1], "block", n_dim_topo);
  block.write(pos.data(), vars.data());
}

void Tecplot_file::write_unstruct(Array<int> elements, Array<double> pos, Array<double> vars)
{
  Unstructured zone(*this, elements, pos, vars);
}

Tecplot_file::~Tecplot_file()
{
  tecFileWriterClose(&file_handle);
}

Tecplot_file::Zone::Zone(Tecplot_file& file, int n_nodes, std::string name_arg)
: file{file}, name{name_arg + std::to_string(file.i_zone++)}, n_nodes{n_nodes}, tecio_zone_index{1}, n_total_vars{file.n_dim + file.n_var},
  var_types(n_total_vars, 2), // declare data type as double
  shared(n_total_vars, 0), // declare no variables shared
  location(n_total_vars, 1), // declare variable location as node-centered
  passive(n_total_vars, 0), // declare no variables passive
  strand_id{file.strand_id++}
{}

Tecplot_file::Zone::~Zone()
{
  tecZoneSetUnsteadyOptions(file.file_handle, tecio_zone_index, file.time, strand_id);
}

void Tecplot_file::Zone::write(const double* pos, const double* vars)
{
  for (int i_dim = 0; i_dim < file.n_dim; ++i_dim) {
    tecZoneVarWriteDoubleValues(file.file_handle, tecio_zone_index, i_dim + 1, 0, n_nodes, pos + i_dim*n_nodes);
  }
  for (int i_var = 0; i_var < file.n_var; ++i_var) {
    tecZoneVarWriteDoubleValues(file.file_handle, tecio_zone_index, file.n_dim + i_var + 1, 0, n_nodes, vars + i_var*n_nodes);
  }
}

Tecplot_file::Structured_block::Structured_block(Tecplot_file& file, int row_size, std::string name_arg, int n_dim_arg)
: Zone{file, math::pow(row_size, n_dim_arg ? n_dim_arg : file.n_dim), name_arg}, n_dim{n_dim_arg ? n_dim_arg : file.n_dim}
{
  tecZoneCreateIJK(file.file_handle, name.c_str(),
                   row_size, (n_dim >= 2) ? row_size : 1, (n_dim >= 3) ? row_size : 1,
                   var_types.data(), shared.data(), location.data(), passive.data(), 0, 0, 0, &tecio_zone_index);
}

Tecplot_file::Line_segments::Line_segments(Tecplot_file& file, int n_segs, int row_size, std::string name_arg)
: Zone{file, row_size*n_segs, name_arg}, n_segs{n_segs}, row_size{row_size}, i_seg{0}
{
  pos_storage.resize(file.n_dim*n_nodes);
  var_storage.resize(file.n_var*n_nodes);

  tecZoneCreateFE(file.file_handle, name.c_str(), us_line_seg, n_nodes, n_segs,
                  var_types.data(), shared.data(), location.data(), passive.data(), 0, 0, 0, &tecio_zone_index);
}

void Tecplot_file::Line_segments::write(const double* pos, const double* vars)
{
  for (int i_node = 0; i_node < row_size; ++i_node) {
    for (int i_dim = 0; i_dim < file.n_dim; ++i_dim) {
      pos_storage[(n_segs*i_dim + i_seg)*row_size + i_node] = pos[i_dim*row_size + i_node];
    }
    for (int i_var = 0; i_var < file.n_var; ++i_var) {
      var_storage[(n_segs*i_var + i_seg)*row_size + i_node] = vars[i_var*row_size + i_node];
    }
  }
  ++i_seg;
}

Tecplot_file::Line_segments::~Line_segments()
{
  Zone::write(pos_storage.data(), var_storage.data());
  std::vector<int> inds;
  for (int i_seg = 0; i_seg < n_segs; ++i_seg) {
    for (int i_elem = 0; i_elem < row_size - 1; ++i_elem) {
      inds.push_back(row_size*i_seg + i_elem);
      inds.push_back(row_size*i_seg + i_elem + 1);
    }
  }
  tecZoneNodeMapWrite32(file.file_handle, tecio_zone_index, 0, 0, inds.size(), inds.data());
}

Tecplot_file::Triangles::Triangles(Tecplot_file& file, int n_triangles, std::string name_arg)
: Zone{file, 3, name_arg}, n_tri{n_triangles}
{
  tecZoneCreateFE(file.file_handle, name.c_str(), us_triangle, 3*n_triangles, n_triangles,
                  var_types.data(), shared.data(), location.data(), passive.data(), 0, 0, 0, &tecio_zone_index);
}

Tecplot_file::Triangles::~Triangles()
{
  std::vector<int> inds;
  for (int i = 0; i < 3*n_tri; ++i) inds.push_back(i);
  tecZoneNodeMapWrite32(file.file_handle, tecio_zone_index, 0, 0, inds.size(), inds.data());
}

Tecplot_file::Unstructured::Unstructured(Tecplot_file& file, Array<Int> elements, Array<double> pos, Array<double> vars, std::string name_arg)
: Zone{file, 0, name_arg}
{
  HEXED_ASSERT(pos.order() == 2, "order of `pos` must be 2");
  HEXED_ASSERT(pos.shape()[0] == file.n_dim, "must have same number of position variables as file dimensionality");
  HEXED_ASSERT(pos(0).same_shape(vars(0)), "`pos` and `vars` must have same number of variables and same order");
  HEXED_ASSERT(vars.shape()[0] == file.n_var, "must have same number of field variables as file");
  const char* name = "unstructured_zone";
  Zone z(file, 3, name);
  int n_vert = pos.shape()[1];
  int n_elem_vert = elements.shape()[1];
  int n_elem = elements.shape()[0];
  int zone_t;
  if      (n_elem_vert == 2 && file.n_dim_topo == 1) zone_t = us_line_seg;
  else if (n_elem_vert == 3 && file.n_dim_topo == 2) zone_t = us_triangle;
  else if (n_elem_vert == 4 && file.n_dim_topo == 2) zone_t = us_quad;
  else if (n_elem_vert == 4 && file.n_dim_topo == 3) zone_t = us_tet;
  else if (n_elem_vert == 8 && file.n_dim_topo == 3) zone_t = us_hex;
  else HEXED_ASSERT(false, "current combination of number of vertices and toplogical dimension did not match any known element type");
  HEXED_ASSERT(zone_t != us_quad && zone_t != us_hex, "vertex permutations for unstructured quad/hex visualization with Tecplot have not been implemented");
  tecZoneCreateFE(file.file_handle, name, zone_t, pos.shape()[1], n_elem,
                  var_types.data(), shared.data(), location.data(), passive.data(), 0, 0, 0, &tecio_zone_index);
  for (int i_dim = 0; i_dim < file.n_dim; ++i_dim) {
    tecZoneVarWriteDoubleValues(file.file_handle, tecio_zone_index, i_dim + 1, 0, n_vert, pos(i_dim).data());
  }
  for (int i_var = 0; i_var < file.n_var; ++i_var) {
    tecZoneVarWriteDoubleValues(file.file_handle, tecio_zone_index, file.n_dim + i_var + 1, 0, n_vert, vars(i_var).data());
  }
  tecZoneNodeMapWrite32(file.file_handle, tecio_zone_index, 0, 0, elements.size(), elements.data());
}

void Tecplot_file::Unstructured::write(const double* pos, const double* vars)
{
  HEXED_ASSERT(false, "don't use `write` for `hexed::Tecplot_file::Unstructured`. Constructor does it for you.");
}

}
#endif
