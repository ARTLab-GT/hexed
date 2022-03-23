#include <Accessible_mesh.hpp>
#include <math.hpp>

namespace cartdg
{

Element_container& Accessible_mesh::container(bool is_deformed)
{
  Element_container* containers [] {&car.elems, &def.elems};
  return *containers[is_deformed];
}

Accessible_mesh::Accessible_mesh(Storage_params params_arg, double root_size) :
  params{params_arg},
  root_sz{root_size},
  car{params, root_sz},
  def{params, root_sz},
  def_as_car{def.elements()},
  elems{car.elements(), def_as_car},
  elem_cons{car.element_connections(),
  def.element_connections()},
  bc_v{bound_conds},
  bound_face_cons{car.bound_face_con_view, def.bound_face_con_view},
  def_face_cons{def.elem_face_con_v, bound_face_cons},
  ref_face_v{car.refined_faces(), def.refined_faces()}
{
  def.face_con_v = def_face_cons;
}

int Accessible_mesh::add_element(int ref_level, bool is_deformed, std::vector<int> position)
{
  return container(is_deformed).emplace(ref_level, position);
}

Element& Accessible_mesh::element(int ref_level, bool is_deformed, int serial_n)
{
  return container(is_deformed).at(ref_level, serial_n);
}

void Accessible_mesh::connect_cartesian(int ref_level, std::array<int, 2> serial_n, Con_dir<Element> direction, std::array<bool, 2> is_deformed)
{
  std::array<Element*, 2> el_ar;
  for (int i_side : {0, 1}) el_ar[i_side] = &element(ref_level, is_deformed[i_side], serial_n[i_side]);
  car.cons.emplace_back(el_ar, direction);
}

void Accessible_mesh::connect_deformed(int ref_level, std::array<int, 2> serial_n, Con_dir<Deformed_element> direction)
{
  if ((direction.i_dim[0] == direction.i_dim[1]) && (direction.face_sign[0] == direction.face_sign[1])) {
    throw std::runtime_error("attempt to connect faces of same sign along same dimension which is forbidden");
  }
  std::array<Deformed_element*, 2> el_ar;
  for (int i_side : {0, 1}) {
    el_ar[i_side] = &def.elems.at(ref_level, serial_n[i_side]);
  }
  def.cons.emplace_back(el_ar, direction);
}

void Accessible_mesh::connect_hanging_cartesian(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial, Con_dir<Element> dir,
                                                bool coarse_face_positive, bool coarse_deformed, bool fine_deformed)
{
  Element* coarse = &element(coarse_ref_level, coarse_deformed, coarse_serial);
  std::vector<Element*> fine;
  for (int fs : fine_serial) fine.push_back(&element(coarse_ref_level + 1, fine_deformed, fs));
  car.ref_face_cons.emplace_back(coarse, fine, dir, !coarse_face_positive);
}

int Accessible_mesh::add_boundary_condition(Boundary_condition* bc)
{
  bound_conds.emplace_back(bc);
  // no reason to delete boundary conditions, so the serial number can just be the index
  return bound_conds.size() - 1;
}

void Accessible_mesh::connect_boundary(int ref_level, bool is_deformed, int element_serial_n, int i_dim, int face_sign, int bc_serial_n)
{
  if (bc_serial_n >= int(bound_conds.size())) throw std::runtime_error("demand for non-existent `Boundary_condition`");
  Boundary_condition& bc {*bound_conds[bc_serial_n]};
  #define EMPLACE(mbt) { \
    mbt.bound_cons.emplace_back(mbt.elems.at(ref_level, element_serial_n), i_dim, face_sign, bc); \
  }
  if (is_deformed) EMPLACE(def)
  else EMPLACE(car)
  #undef EMPLACE
}

Mesh::Connection_validity Accessible_mesh::valid()
{
  // maps a pointer to an element face to the number of connections pointing to that face
  std::map<double*, int> n_connections;
  // add an entry for each face and initialize it to 0
  auto& elems = elements();
  int n_dof_face = params.n_dof()/params.row_size;
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* faces = elems[i_elem].face();
    for (int i_face = 0; i_face < params.n_dim*2; ++i_face) {
      n_connections[faces + i_face*n_dof_face] = 0;
    }
  }
  // count up the number of connections for each face
  #define COUNT_CONS(mbt) \
  { \
    auto& cons = mbt.face_connections(); \
    for (int i_con = 0; i_con < cons.size(); ++i_con) { \
      for (int i_side = 0; i_side < 2; ++i_side) { \
        double* face = cons[i_con].face(i_side); \
        if (n_connections.count(face)) { \
          ++n_connections.at(face); \
        } \
      } \
    } \
  }
  COUNT_CONS(car)
  COUNT_CONS(def)
  auto& ref_faces = refined_faces();
  for (int i_ref = 0; i_ref < ref_faces.size(); ++i_ref) {
    ++n_connections.at(ref_faces[i_ref].coarse_face());
  }
  // count up the number of faces with problems
  int n_missing = 0;
  int n_duplicates = 0;
  for (auto pair : n_connections) {
    if (pair.second == 0) ++n_missing;
    n_duplicates += std::max(pair.second - 1, 0);
  }
  return {n_duplicates, n_missing};
}

}
