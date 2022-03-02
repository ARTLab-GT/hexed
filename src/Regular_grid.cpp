#include <Regular_grid.hpp>
#include <get_mcs_convective.hpp>
#include <get_write_face.hpp>
#include <get_prolong.hpp>
#include <get_neighbor_convective.hpp>
#include <get_restrict.hpp>
#include <get_local_convective.hpp>
#include <get_gbc_convective.hpp>
#include <get_req_visc_regular_convective.hpp>
#include <get_write_face_scalar.hpp>
#include <get_neighbor_gradient.hpp>
#include <get_local_gradient.hpp>
#include <get_write_face_n_dim.hpp>
#include <get_neighbor_av.hpp>
#include <get_local_av.hpp>
#include <math.hpp>

namespace cartdg
{

Regular_grid::Regular_grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: Grid (n_var_arg, n_dim_arg, n_elem_arg, mesh_size_arg, basis_arg), elements{}
{
  for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
    elements.emplace_back(new Element {storage_params});
  }
  for (int i_dim = 0; i_dim < 3*n_dim; ++i_dim) {
    elem_cons.push_back({});
    ref_faces.push_back({});
  }
}

Element& Regular_grid::element(int i_elem)
{
  return *elements[i_elem];
}

elem_con Regular_grid::connection(int i_dim, int i_con)
{
  return elem_cons[i_dim][i_con];
}

int Regular_grid::n_con(int i_dim)
{
  return int(elem_cons[i_dim].size());
}

std::vector<double> Regular_grid::get_pos(int i_elem)
{
  std::vector<double> elem_pos (n_qpoint*n_dim, 0.);
  std::vector<int> indices;
  populate_slice(elem_pos, indices, i_elem);
  return elem_pos;
}

double Regular_grid::stable_time_step(double cfl_by_stable_cfl, Kernel_settings& settings)
{
  double cfl = cfl_by_stable_cfl*get_stable_cfl();
  settings.i_read = i_read;
  return cfl*mesh_size/get_mcs_convective(n_dim, basis.row_size)(elements, settings);
}

void Regular_grid::execute_write_face(Kernel_settings& settings)
{
  settings.i_read = i_read;
  get_write_face(n_dim, basis.row_size)(elements, basis, settings);
  // it's important that the `write_face` of lower-ref-level grids has already happened
  // if no refined faces, don't even call the function (to allow to use basis without prolong method as long as grid doesn't have hanging nodes)
  if (ref_faces[0].size()) get_prolong(n_dim, basis.row_size)(ref_faces, basis, settings);
}

void Regular_grid::execute_neighbor(Kernel_settings& settings)
{
  get_neighbor_convective(n_dim, basis.row_size)(elem_cons, settings);
  get_gbc_convective(n_dim, basis.row_size)(element_gbcs, basis, settings);
  if (ref_faces[0].size()) get_restrict(n_dim, basis.row_size)(ref_faces, basis, settings);
}

void Regular_grid::execute_local(Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_local_convective(n_dim, basis.row_size)(elements, basis, settings);
}

void Regular_grid::execute_req_visc(Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_req_visc_regular_convective(n_dim, basis.row_size)(elements, basis, settings);
}

void Regular_grid::execute_write_face_gradient(int i_var, Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_write_face_scalar(n_dim, basis.row_size)(elements, i_var, basis, settings);
}

void Regular_grid::execute_neighbor_gradient(int i_var, Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_neighbor_gradient(n_dim, basis.row_size)(elem_cons, i_var, settings);
}

void Regular_grid::execute_local_gradient(int i_var, Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_local_gradient(n_dim, basis.row_size)(elements, i_var, basis, settings);
}

void Regular_grid::execute_write_face_av(int i_var, Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_write_face_n_dim(n_dim, basis.row_size)(elements, basis, settings);
}

void Regular_grid::execute_neighbor_av(int i_var, Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_neighbor_av(n_dim, basis.row_size)(elem_cons, i_var, settings);
}

void Regular_grid::execute_local_av(int i_var, Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_local_av(n_dim, basis.row_size)(elements, i_var, basis, settings);
}

int Regular_grid::add_element(std::vector<int> position)
{
  elements.emplace_back(new Element {storage_params, position, mesh_size});
  return Grid::add_element(position);
}

void Regular_grid::add_element_gbc(int i_elem, Ghost_boundary_condition& gbc)
{
  element_gbcs.emplace_back(element(i_elem), gbc);
}

void Regular_grid::add_connection(int i_elem0, int i_elem1, int i_dim)
{
  add_connection(elements[i_elem0].get(), elements[i_elem1].get(), i_dim);
}

void Regular_grid::add_connection(Element* elem0, Element* elem1, int i_dim)
{
  const int fs {n_var*n_qpoint/basis.row_size};
  elem_con con {elem0->face() + (2*i_dim + 1)*fs, elem1->face() + 2*i_dim*fs};
  elem_cons[i_dim].push_back(con);
  int stride = custom_math::pow(2, n_dim - 1 - i_dim);
  for (int i_vert = 0; i_vert < n_vertices/2; ++i_vert) {
    int i_col = i_vert/stride*stride*2 + i_vert%stride;
    elem0->vertex(i_col + stride).eat(elem1->vertex(i_col));
  }
}

void Regular_grid::connect_refined(Element* coarse, std::vector<Element*> fine, int i_dim, bool is_positive)
{
  // create a `Hanging_node_matcher` to manage shareable data at the hanging vertices
  hanging_matchers.emplace_back(fine, i_dim, !is_positive); // `is_positive` refers to coarse not fine
  const int fs {n_var*n_qpoint/basis.row_size};
  if (fine.size() != unsigned(n_vertices)/2) throw std::runtime_error("Wrong number of fine elements in hanging-node connection.");
  // create a `Refined_face` to store the mortar faces
  ref_faces[i_dim].emplace_back(new Refined_face {storage_params, coarse->face() + (2*i_dim + is_positive)*fs});
  int vertex_stride = custom_math::pow(2, n_dim - 1 - i_dim);
  std::vector<Element*> face_elements;
  for (int i_face = 0; i_face < n_vertices/2; ++i_face) {
    // connect the mortar faces to the fine faces
    elem_con con {ref_faces[i_dim].back()->fine_face(i_face), fine[i_face]->face() + (2*i_dim + 1 - is_positive)*fs};
    if (!is_positive) std::swap(con[0], con[1]);
    elem_cons[i_dim].push_back(con);
    // combine coarse vertices with matching fine vertices
    int vertex_col = i_face/vertex_stride*vertex_stride*2 + i_face%vertex_stride;
    coarse->vertex(vertex_col + is_positive*vertex_stride).eat(fine[i_face]->vertex(vertex_col + (1 - is_positive)*vertex_stride));
  }
}

void Regular_grid::auto_connect(std::vector<int> periods)
{
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int j_elem = i_elem + 1; j_elem < n_elem; ++j_elem)
    {
      int pos_diff [n_dim];
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        pos_diff[i_dim] = pos[n_dim*j_elem + i_dim] - pos[n_dim*i_elem + i_dim];
        // assumes period[i_dim] >= greatest distance between elements in this dimension
        if (periods[i_dim] > 0)
        {
          pos_diff[i_dim] = (pos_diff[i_dim] + periods[i_dim] + 1)%periods[i_dim] - 1;
        }
      }
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        bool is_same_row = true;
        for (int j_dim = 0; j_dim < n_dim; ++j_dim)
        {
          if ((j_dim != i_dim) && (pos_diff[j_dim] != 0)) is_same_row = false;
        }
        if (is_same_row)
        {
          if      (pos_diff[i_dim] ==  1) add_connection(i_elem, j_elem, i_dim);
          else if (pos_diff[i_dim] == -1) add_connection(j_elem, i_elem, i_dim);
        }
      }
    }
  }
}

void Regular_grid::auto_connect()
{
  std::vector<int> periods;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) periods.push_back(0);
  auto_connect(periods);
}

void Regular_grid::populate_slice(std::vector<double>& elem_pos, std::vector<int> indices, int i_elem)
{
  if ((int)indices.size() < n_dim)
  {
    indices.push_back(0);
    for (int i = 0; i < basis.row_size; ++i) {
      indices.back() = i;
      populate_slice(elem_pos, indices, i_elem);
    }
  }
  else
  {
    int i_flat = 0;
    int stride = n_qpoint;
    for (auto i : indices) {
      stride /= basis.row_size;
      i_flat += i*stride;
    }
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      elem_pos[i_flat + n_qpoint*i_dim] = (basis.node(indices[i_dim])
                                           + pos[i_elem*n_dim + i_dim])*mesh_size;
    }
  }
}

}
