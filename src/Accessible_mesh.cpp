#include <Accessible_mesh.hpp>
#include <math.hpp>
#include <erase_if.hpp>

namespace cartdg
{

Element_container& Accessible_mesh::container(bool is_deformed)
{
  Element_container* containers [] {&car.elems, &def.elems};
  return *containers[is_deformed];
}

Accessible_mesh::Accessible_mesh(Storage_params params_arg, double root_size) :
  params{params_arg},
  n_vert{custom_math::pow(2, params.n_dim)},
  root_sz{root_size},
  car{params, root_sz},
  def{params, root_sz},
  def_as_car{def.elements()},
  elems{car.elements(), def_as_car},
  elem_cons{car.element_connections(), def.element_connections()},
  bc_v{bound_conds},
  bound_face_cons{car.bound_face_con_view, def.bound_face_con_view},
  bound_cons{car.boundary_connections(), def.boundary_connections()},
  def_face_cons{def.elem_face_con_v, bound_face_cons},
  ref_face_v{car.refined_faces(), def.refined_faces()},
  matcher_v{car.hanging_vertex_matchers(), def.hanging_vertex_matchers()}
{
  def.face_con_v = def_face_cons;
}

int Accessible_mesh::add_element(int ref_level, bool is_deformed, std::vector<int> position)
{
  int sn = container(is_deformed).emplace(ref_level, position);
  Element& elem = element(ref_level, is_deformed, sn);
  for (int i_vert = 0; i_vert < n_vert; ++i_vert) vert_ptrs.emplace_back(elem.vertex(i_vert));
  return sn;
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
  car.ref_face_cons.emplace_back(new Refined_connection<Element> {coarse, fine, dir, !coarse_face_positive});
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
  auto& elems = elements();
  const int n_faces = 2*params.n_dim;
  // initialize number of connections of each face to 0
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    for (int i_face = 0; i_face < n_faces; ++i_face) {
      elems[i_elem].face_record[i_face] = 0;
    }
  }
  // count up the number of connections for each face
  car.record_connections();
  def.record_connections();
  // count up the number of faces with problems
  int n_missing = 0;
  int n_duplicates = 0;
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    for (int i_face = 0; i_face < n_faces; ++i_face) {
      int rec = elems[i_elem].face_record[i_face];
      if (rec == 0) ++n_missing;
      if (rec > 1) n_duplicates += rec - 1;
    }
  }
  return {n_duplicates, n_missing};
}

Accessible_mesh::vertex_view Accessible_mesh::vertices()
{
  erase_if(vert_ptrs, &Vertex::Non_transferable_ptr::is_null);
  return vert_ptrs;
}

// helper objects for `extrude`
struct Empty_face
{
  Deformed_element& elem;
  int i_dim;
  int face_sign;
};
struct Connection_plan
{
  int ref_level;
  std::array<int, 2> serial_ns;
  Con_dir<Deformed_element> dir;
};
bool aligned_same_dim(Con_dir<Deformed_element> dir, std::array<int, 2> extrude_dim)
{
  return (dir.i_dim[0] == dir.i_dim[1]) && (dir.face_sign[0] != dir.face_sign[1])
         && (extrude_dim[0] == extrude_dim[1]);
}
bool aligned_different_dim(Con_dir<Deformed_element> dir, std::array<int, 2> extrude_dim)
{
  return (dir.i_dim[0] == extrude_dim[1]) && (dir.i_dim[1] == extrude_dim[0]);
}

void Accessible_mesh::extrude()
{
  const int n_faces = 2*params.n_dim;
  // initialize number of connections of each face to 0
  {
    auto& elems = elements();
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      for (int i_face = 0; i_face < n_faces; ++i_face) {
        elems[i_elem].face_record[i_face] = 0;
      }
    }
  }
  // initialize vertex records to empty
  {
    auto verts = vertices();
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
      verts[i_vert].record.clear();
    }
  }
  // count up the number of connections for each face
  car.record_connections();
  def.record_connections();
  // record which faces have boundary conditions
  auto& bc_cons {boundary_connections()};
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    auto& con {bc_cons[i_con]};
    // find the index of the boundary condition
    int i_bc = 0;
    for (; &boundary_conditions()[i_bc] != con.boundary_condition(); ++i_bc);
    // no face has more than one connections (I really hope!) so use numbers greater than one to indentify boundary conditions
    con.element().face_record[2*con.i_dim() + con.inside_face_sign()] = 2 + i_bc;
  }
  // decide which faces to extrude from
  std::vector<Empty_face> empty_faces;
  auto& elems = def.elements();
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
      for (int face_sign = 0; face_sign < 2; ++face_sign) {
        const int i_face = 2*i_dim + face_sign;
        auto& elem {elems[i_elem]};
        if (elem.face_record[i_face] == 0) {
          empty_faces.push_back({elem, i_dim, face_sign});
        }
      }
    }
  }
  // create extruded elements
  const int n_record = 4;
  for (auto face : empty_faces)
  {
    auto nom_pos = face.elem.nominal_position();
    nom_pos[face.i_dim] += 2*face.face_sign - 1;
    const int ref_level = face.elem.refinement_level();
    int sn = add_element(ref_level, true, nom_pos);
    Con_dir<Deformed_element> dir {{face.i_dim, face.i_dim}, {!face.face_sign, bool(face.face_sign)}};
    auto& elem = def.elems.at(ref_level, sn);
    std::array<Deformed_element*, 2> el_arr {&elem, &face.elem};
    def.cons.emplace_back(el_arr, dir);
    // record the faces that still need to be connected at a vertex which is guaranteed to be shared with prospective neighbors
    for (int j_dim = face.i_dim + 1; j_dim%params.n_dim != face.i_dim; ++j_dim)
    {
      j_dim = j_dim%params.n_dim;
      for (int face_sign = 0; face_sign < 2; ++face_sign)
      {
        int face_rec = face.elem.face_record[2*j_dim + face_sign];
        if (face_rec >= 2)
        {
          def.bound_cons.emplace_back(elem, j_dim, face_sign, boundary_conditions()[face_rec - 2]);
        }
        else
        {
          // record data at vertex which is on the face to be connected, on the face which was extruded from,
          // and if applicable has the minimum index to satisfy the above consitions.
          int i_vert =   face_sign*custom_math::pow(2, params.n_dim - 1 - j_dim)
                       + (1 - face.face_sign)*custom_math::pow(2, params.n_dim - 1 - face.i_dim);
          auto& record = elem.vertex(i_vert).record;
          // which element it is
          record.push_back(ref_level);
          record.push_back(sn);
          // which face needs to be connected
          record.push_back(2*j_dim + face_sign);
          // extrusion direction for deciding which face connections are valid
          record.push_back(face.i_dim);
        }
      }
    }
  }
  // plan connections to make (don't make them yet, because that could result in `eat`ing vertices which have not been visited,
  // and ultimately dereferencing null pointers)
  std::vector<Connection_plan> con_plans;
  {
    auto verts = vertices(); // note: need to rebuild vertex vector because face connections above have `eat`en vertices
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert)
    {
      auto& vert {verts[i_vert]};
      // first make connections where dimension matches and then make connections among differing dimensions.
      // This order prevents incorrect connections where both same-dimension and different-dimension candidates are available.
      for (bool (*aligned)(Con_dir<Deformed_element>, std::array<int, 2>) : {&aligned_same_dim, &aligned_different_dim})
      {
        // iterate through every possible pair of records created by an extruded elements above
        for (int i_record = 0; i_record < int(vert.record.size()); i_record += n_record)
        {
          for (unsigned j_record = i_record + n_record; j_record < vert.record.size(); j_record += n_record)
          {
            int ref_level = vert.record[i_record];
            std::array<int, 2> sn {vert.record[i_record + 1], vert.record[j_record + 1]};
            Con_dir<Deformed_element> dir {{     vert.record[i_record + 2]/2 ,      vert.record[j_record + 2]/2},
                                           {bool(vert.record[i_record + 2]%2), bool(vert.record[j_record + 2]%2)}};
            // only connect elements that are suitably positioned.
            // This prevents incorrect connections at places like a 3D corner where there are many (incorrect) candidates available
            if (aligned(dir, {vert.record[i_record + 3], vert.record[j_record + 3]})) {
              if (vert.record[j_record] != ref_level) throw std::runtime_error("attempt to connect extruded elements of different ref level");
              con_plans.push_back({ref_level, sn, dir}); // add prospective connection
              // erase unconnected face record to prevent duplicate connections
              vert.record.erase(vert.record.begin() + j_record, vert.record.begin() + j_record + n_record);
              vert.record.erase(vert.record.begin() + i_record, vert.record.begin() + i_record + n_record);
              i_record -= n_record; // move index to account for erased elements
              break; // since we found a match, we can move on to the next `i_record`
            }
          }
        }
      }
    }
  }
  // create the planned connections
  for (auto con_plan : con_plans) {
    connect_deformed(con_plan.ref_level, con_plan.serial_ns, con_plan.dir);
  }
}

void Accessible_mesh::connect_rest(int bc_sn)
{
  auto& bc {boundary_condition(bc_sn)};
  car.connect_rest(bc);
  def.connect_rest(bc);
}

}
