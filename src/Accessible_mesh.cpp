#include <Accessible_mesh.hpp>
#include <math.hpp>
#include <Row_index.hpp>
#include <erase_if.hpp>
#include <utils.hpp>

namespace hexed
{

Element_container& Accessible_mesh::container(bool is_deformed)
{
  Element_container* containers [] {&car.elems, &def.elems};
  return *containers[is_deformed];
}

template<> Mesh_by_type<         Element>& Accessible_mesh::mbt() {return car;}
template<> Mesh_by_type<Deformed_element>& Accessible_mesh::mbt() {return def;}

Accessible_mesh::Accessible_mesh(Storage_params params_arg, double root_size_arg) :
  params{params_arg},
  n_vert{math::pow(2, params.n_dim)},
  root_sz{root_size_arg},
  car{params, root_sz},
  def{params, root_sz},
  def_as_car{def.elements()},
  elems{car.elements(), def_as_car},
  elem_cons{car.element_connections(), def.element_connections()},
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

void Accessible_mesh::connect_hanging(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial, Con_dir<Deformed_element> dir,
                                      bool coarse_deformed, std::vector<bool> fine_deformed, std::array<bool, 2> stretch)
{
  bool is_car = !coarse_deformed;
  for (bool fine_def : fine_deformed) is_car = (is_car||!fine_def);
  if (is_car)
  {
    Element* coarse = &element(coarse_ref_level, coarse_deformed, coarse_serial);
    std::vector<Element*> fine;
    for (int i_fine = 0; i_fine < n_vert/2; ++i_fine) {
      fine.push_back(&element(coarse_ref_level + 1, fine_deformed[i_fine], fine_serial[i_fine]));
    }
    if ((dir.i_dim[0] != dir.i_dim[1]) || (dir.face_sign[0] == dir.face_sign[1])) {
      throw std::runtime_error("attempted to form a cartesian hanging-node connection with incompatible `Con_dir`.");
    }
    car.ref_face_cons[params.n_dim - 1].emplace_back(new Refined_connection<Element> {coarse, fine, {dir.i_dim[0]}, dir.face_sign[1]});
  }
  else
  {
    Deformed_element* coarse = &def.elems.at(coarse_ref_level, coarse_serial);
    std::vector<Deformed_element*> fine;
    for (unsigned i_fine = 0; i_fine < fine_serial.size(); ++i_fine) {
      fine.push_back(&def.elems.at(coarse_ref_level + 1, fine_serial[i_fine]));
    }
    def.ref_face_cons[math::log(2, fine_serial.size())].emplace_back(new Refined_connection<Deformed_element> {coarse, fine, dir, false, stretch});
  }
}

int Accessible_mesh::add_boundary_condition(Flow_bc* flow_bc, Mesh_bc* mesh_bc)
{
  bound_conds.push_back({std::unique_ptr<Flow_bc>{flow_bc}, std::unique_ptr<Mesh_bc>{mesh_bc}});
  // no reason to delete boundary conditions, so the serial number can just be the index
  return bound_conds.size() - 1;
}

void Accessible_mesh::connect_boundary(int ref_level, bool is_deformed, int element_serial_n, int i_dim, int face_sign, int bc_serial_n)
{
  if (bc_serial_n >= int(bound_conds.size())) throw std::runtime_error("demand for non-existent `Boundary_condition`");
  #define EMPLACE(mbt) { \
    mbt.bound_cons.emplace_back(mbt.elems.at(ref_level, element_serial_n), i_dim, face_sign, bc_serial_n); \
  }
  if (is_deformed) EMPLACE(def)
  else EMPLACE(car)
  #undef EMPLACE
}

void Accessible_mesh::disconnect_boundary(int bc_sn)
{
  auto is_doomed = [bc_sn](Boundary_connection& con){return con.bound_cond_serial_n() == bc_sn;};
  erase_if(car.bound_cons, is_doomed);
  erase_if(def.bound_cons, is_doomed);
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

//! \cond helper classes and functions for Accessible_mesh::extrude
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
struct Refined_connection_plan
{
  int coarse_ref;
  int coarse_sn;
  std::vector<int> fine_sn;
  Con_dir<Deformed_element> dir;
  std::array<bool, 2> stretch;
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
//! \endcond

void Accessible_mesh::extrude(bool collapse, double offset)
{
  const int nd = params.n_dim;
  const int n_faces = 2*nd;
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
    // no face has more than one connections (I really hope!) so use numbers greater than one to indentify boundary conditions
    con.element().face_record[2*con.i_dim() + con.inside_face_sign()] = 2 + con.bound_cond_serial_n();
  }
  // decide which faces to extrude from
  std::vector<Empty_face> empty_faces;
  auto& elems = def.elements();
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
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
    if (collapse) {
      int stride = math::pow(2, nd - 1 - face.i_dim);
      for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
        int i_collapse = i_vert + (face.face_sign - (i_vert/stride)%2)*stride;
        elem.vertex(i_vert).pos = face.elem.vertex(i_collapse).pos;
      }
    }
    std::array<Deformed_element*, 2> el_arr {&elem, &face.elem};
    extrude_cons.push_back(def.cons.size());
    def.cons.emplace_back(el_arr, dir);
    // record the faces that still need to be connected at a vertex which is guaranteed to be shared with prospective neighbors
    for (int j_dim = face.i_dim + 1; j_dim%nd != face.i_dim; ++j_dim)
    {
      j_dim = j_dim%nd;
      for (int face_sign = 0; face_sign < 2; ++face_sign)
      {
        int face_rec = face.elem.face_record[2*j_dim + face_sign];
        if (face_rec >= 2)
        {
          // if parent element has boundary connections on other faces
          def.bound_cons.emplace_back(elem, j_dim, face_sign, face_rec - 2);
        }
        else
        {
          // record data at vertex which is on the face to be connected, on the face which was extruded from,
          // and if applicable has the minimum index to satisfy the above consitions.
          int i_vert =   face_sign*math::pow(2, nd - 1 - j_dim)
                       + (1 - face.face_sign)*math::pow(2, nd - 1 - face.i_dim);
          auto& record = elem.vertex(i_vert).record;
          // which element it is
          record.push_back(ref_level);
          record.push_back(sn);
          // which face needs to be connected
          record.push_back(2*j_dim + face_sign);
          // extrusion direction for deciding which face connections are valid
          record.push_back(2*face.i_dim + face.face_sign);
        }
      }
    }
    // readjust face node adjustments to account for offset
    if (offset > 0) {
      double* node_adj [] {face.elem.node_adjustments(), elem.node_adjustments()};
      int nfq = params.n_qpoint()/params.row_size;
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
        double* na [2][2];
        for (int i = 0; i < 2; ++i) {
          for (int face_sign = 0; face_sign < 2; ++face_sign) {
            na[i][face_sign] = node_adj[i] + (2*face.i_dim + (face_sign == face.face_sign))*nfq + i_qpoint;
          }
        }
        *na[1][1] = *na[0][1];
        *na[0][1] = *na[1][0] = offset**na[0][0] + (1 - offset)**na[0][1];
        for (int sign = 0; sign < 2; ++sign) {
          *na[0][sign] /= 1 - offset;
          *na[1][sign] /= offset;
        }
      }
    }
  }

  // perform offset
  if (offset > 0) {
    // loop through all vertices on faces that used to be exposed
    for (auto face : empty_faces) {
      for (Row_index ind(nd, 2, face.i_dim); ind; ++ind) {
        Vertex& vert = face.elem.vertex(ind.i_qpoint(face.face_sign)); // vertex that has to be moved
        Vertex& off_vert = face.elem.vertex(ind.i_qpoint(!face.face_sign)); // vertex to move toward
        // use an extra zero appended to the vertex record to mark vertices that have already been moved,
        // so that neighboring elements don't move the same vertices repeatedly
        if (vert.record.size()%4 == 0) {
          vert.record.push_back(0);
          // move vertex
          for (int i_dim = 0; i_dim < nd; ++i_dim) {
            vert.pos[i_dim] = (1 - offset)*vert.pos[i_dim] + offset*off_vert.pos[i_dim];
          }
        }
      }
    }
    // delete the extra 0s we added to the vertex record so as not to interfere with the connection process below
    for (auto face : empty_faces) {
      for (Row_index ind(nd, 2, face.i_dim); ind; ++ind) {
        Vertex& vert = face.elem.vertex(ind.i_qpoint(face.face_sign));
        if (vert.record.size()%4 == 1) vert.record.pop_back();
      }
    }
  }

  // plan connections to make (don't make them yet, because that could result in `eat`ing vertices which have not been visited,
  // and ultimately dereferencing null pointers)
  std::vector<Connection_plan> con_plans;
  std::vector<Refined_connection_plan> ref_con_plans;

  #define VERTEX_LOOP(code) \
    /* connections without hanging nodes */ \
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert) \
    { \
      auto& vert {verts[i_vert]}; \
      /* first make connections where dimension matches and then make connections among differing dimensions. */ \
      /* This order prevents incorrect connections where both same-dimension and different-dimension candidates are available. */ \
      for (bool (*aligned)(Con_dir<Deformed_element>, std::array<int, 2>) : {&aligned_same_dim, &aligned_different_dim}) \
      { \
        /* iterate through every possible pair of records created by an extruded elements above */ \
        for (int i_record = 0; i_record < int(vert.record.size()); i_record += n_record) \
        { \
          for (int j_record = i_record + n_record; j_record < int(vert.record.size()); j_record += n_record) \
          { \
            int ref_level = vert.record[i_record]; \
            Con_dir<Deformed_element> dir {{     vert.record[i_record + 2]/2 ,      vert.record[j_record + 2]/2}, \
                                           {bool(vert.record[i_record + 2]%2), bool(vert.record[j_record + 2]%2)}}; \
            /* only connect elements that are suitably positioned. */ \
            /* This prevents incorrect connections at places like a 3D corner where there are many (incorrect) candidates available */ \
            if (aligned(dir, {vert.record[i_record + 3]/2, vert.record[j_record + 3]/2})) { \
              code \
            } \
          } \
        } \
      } \
    } \

  #define CONNECT_SAME_CONFORMING \
    if (vert.record[j_record] == ref_level) { \
      std::array<int, 2> sn {vert.record[i_record + 1], vert.record[j_record + 1]}; \
      con_plans.push_back({ref_level, sn, dir}); /* add prospective connection */ \
      /* erase unconnected face record to prevent duplicate connections */ \
      vert.record.erase(vert.record.begin() + j_record, vert.record.begin() + j_record + n_record); \
      vert.record.erase(vert.record.begin() + i_record, vert.record.begin() + i_record + n_record); \
      i_record -= n_record; /* move index to account for erased elements */ \
      break; /* since we found a match, we can move on to the next `i_record` */ \
    } \

  #define CONNECT_HANGING \
    if (vert.record[j_record] != ref_level) { \
      if (std::abs(vert.record[j_record] - ref_level) != 1) { \
        throw std::runtime_error("ref levels of neighboring extrusion faces differ by > 1"); \
      } \
      int which_fine = vert.record[j_record] > vert.record[i_record]; \
      int rec [] {i_record, j_record}; \
      /* ref_level, sn, i_face, extrude_dim */ \
      std::vector<int> sn(n_vert/4); \
      int* vert_rec = vert.record.data() + rec[which_fine]; \
      sn[0] = vert_rec[1]; \
      int stretch_dim = 0; \
      int face_dim = vert_rec[2]/2; \
      int face_sign = vert_rec[2]%2; \
      int extr_dim = vert_rec[3]/2; \
      int extr_sign = vert_rec[3]%2; \
      if (nd == 3) { \
        auto& elem = def.elems.at(vert_rec[0], sn[0]); \
        int free_dim = 3 - face_dim - extr_dim; \
        int iv = face_sign*math::pow(2, nd - 1 - face_dim) \
                 + (1 - extr_sign)*math::pow(2, nd - 1 - extr_dim) \
                 + math::pow(2, nd - 1 - free_dim); \
        Vertex& fine_vert = elem.vertex(iv); \
        /* vertex should have exactly one record */ \
        sn[1] = fine_vert.record[1]; \
        HEXED_ASSERT(fine_vert.record.size() == n_record, \
                     format_str(100, "fine vertex has %li != 1 records (position = [%e, %e, %e])", \
                                fine_vert.record.size()/n_record, \
                                fine_vert.pos[0], fine_vert.pos[1], fine_vert.pos[2]).c_str()); \
        stretch_dim = extr_dim > free_dim; \
        fine_vert.record.erase(fine_vert.record.begin(), fine_vert.record.begin() + n_record); \
      } \
      std::array<bool, 2> stretch {false, false}; \
      stretch[stretch_dim] = true; \
      int* coarse_rec = vert.record.data() + rec[!which_fine]; \
      ref_con_plans.push_back({coarse_rec[0], coarse_rec[1], sn, {{coarse_rec[2]/2, face_dim}, {bool(coarse_rec[2]%2), bool(face_sign)}}, stretch}); \
      /* erase unconnected face record to prevent duplicate connections */ \
      vert.record.erase(vert.record.begin() + j_record, vert.record.begin() + j_record + n_record); \
      vert.record.erase(vert.record.begin() + i_record, vert.record.begin() + i_record + n_record); \
      i_record -= n_record; /* move index to account for erased elements */ \
      break; /* since we found a match, we can move on to the next `i_record` */ \
    } \

  {
    auto verts = vertices();
    VERTEX_LOOP(CONNECT_SAME_CONFORMING)
    VERTEX_LOOP(CONNECT_HANGING)
  }
  // create the planned connections
  for (auto con_plan : con_plans) {
    connect_deformed(con_plan.ref_level, con_plan.serial_ns, con_plan.dir);
  }
  for (auto ref_plan : ref_con_plans) {
    connect_hanging(ref_plan.coarse_ref, ref_plan.coarse_sn, ref_plan.fine_sn, ref_plan.dir, true, std::vector<bool>(n_vert/4, true), ref_plan.stretch);
  }
  con_plans.clear();
  {
    auto verts = vertices(); // note: need to rebuild vertex vector because face connections above have `eat`en vertices
    VERTEX_LOOP(CONNECT_SAME_CONFORMING)
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
      HEXED_ASSERT(verts[i_vert].record.empty(), format_str(100, "vertex detected with %lu unprocessed connection requests",
                                                            verts[i_vert].record.size()/n_record).c_str());
    }
  }
  for (auto con_plan : con_plans) {
    connect_deformed(con_plan.ref_level, con_plan.serial_ns, con_plan.dir);
  }
}

void Accessible_mesh::connect_rest(int bc_sn)
{
  auto& elem_seq = elements();
  // locate unconnected faces
  for (int i_elem = 0; i_elem < elem_seq.size(); ++i_elem) {
    for (int i_face = 0; i_face < 2*params.n_dim; ++i_face) {
      elem_seq[i_elem].face_record[i_face] = 0;
    }
  }
  car.record_connections();
  def.record_connections();
  car.connect_empty(bc_sn);
  def.connect_empty(bc_sn);
}

std::vector<Mesh::elem_handle> Accessible_mesh::elem_handles()
{
  std::vector<Mesh::elem_handle> handles;
  for (bool is_deformed : {0, 1}) {
    for (auto ref_sn : container(is_deformed).elem_handles()) {
      handles.push_back({ref_sn[0], is_deformed, ref_sn[1]});
    }
  }
  return handles;
}

Element& Accessible_mesh::add_elem(bool is_deformed, Tree& t)
{
  auto np = t.coordinates();
  int sn = add_element(t.refinement_level(), is_deformed, std::vector<int>(np.begin(), np.end()));
  auto& elem = element(t.refinement_level(), is_deformed, sn);
  elem.record = sn; // put the serial number in the record so it can be used for connections
  elem.tree = &t;
  t.elem = &elem;
  return elem;
}

void Accessible_mesh::add_tree(std::vector<int> serial_numbers)
{
  HEXED_ASSERT(int(serial_numbers.size()) == 2*params.n_dim, "`serial_numbers` has wrong number of elements");
  HEXED_ASSERT(!tree, "each `Mesh` may only contain one tree");
  tree.reset(new Tree(params.n_dim, root_sz));
  tree_bcs = serial_numbers;
  auto& elem = add_elem(false, *tree);
  int sn = elem.record;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      connect_boundary(0, 0, sn, i_dim, sign, serial_numbers[2*i_dim + sign]);
    }
  }
}

template<typename element_t> void Accessible_mesh::connect_new(int start_at)
{
  auto& m = mbt<element_t>();
  auto elems = m.elems.elements();
  int nd = params.n_dim;
  auto connect_refined = [&](Element& elem, int i_dim, int sign, std::vector<Tree*> neighbors) {
    HEXED_ASSERT(int(neighbors.size()) == math::pow(2, nd - 1), format_str(100, "bad number of neighbors %lu (thanks for nothing, ref level smoother)", neighbors.size()))
    bool is_def = true;
    for (Tree* neighbor : neighbors) {
      HEXED_ASSERT(neighbor->elem, "hanging-node connection with nonexistant elements");
      is_def = is_def && neighbor->elem->get_is_deformed();
    }
    if (is_def) {
    } else {
      std::vector<Element*> fine;
      for (Tree* neighbor : neighbors) fine.push_back(neighbor->elem);
      car.ref_face_cons[nd - 1].emplace_back(new Refined_connection<Element>(&elem, fine, {i_dim}, sign));
    }
  };
  for (int i_elem = start_at; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    if (elem.tree) {
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        for (int sign = 0; sign < 2; ++sign) {
          if (!elem.faces[2*i_dim + sign]) { // if this face hasn't been connected already
            Eigen::VectorXi direction = Eigen::VectorXi::Zero(nd);
            direction(i_dim) = math::sign(sign);
            auto neighbors = elem.tree->find_neighbors(direction);
            if (neighbors.empty()) m.bound_cons.emplace_back(elem, i_dim, sign, tree_bcs[2*i_dim + sign]);
            else if (neighbors.size() == 1) {
              if (neighbors[0]->elem) {
                auto& other = *neighbors[0]->elem;
                if (other.refinement_level() == elem.refinement_level()) {
                  if (elem.get_is_deformed() && other.get_is_deformed()) ;
                  else {
                    std::array<Element*, 2> el_ar;
                    el_ar[!sign] = &elem;
                    el_ar[sign] = &other;
                    car.cons.emplace_back(el_ar, Con_dir<Element>{i_dim});
                  }
                } else {
                  bool is_min_corner = true;
                  for (int j_dim = 0; j_dim < nd; ++j_dim) if (j_dim != i_dim) {
                    is_min_corner = is_min_corner && elem.tree->coordinates()[j_dim]%2 == 0;
                  }
                  if (is_min_corner) connect_refined(other, i_dim, sign, other.tree->find_neighbors(-direction));
                }
              }
            } else {
              connect_refined(elem, i_dim, !sign, neighbors);
            }
          }
        }
      }
    }
  }
}

void Accessible_mesh::refine_by_record(bool is_deformed, int start, int end)
{
  auto& elems = container(is_deformed).element_view();
  for (int i_elem = start; i_elem < end; ++i_elem) {
    auto& elem = elems[i_elem];
    if (elem.tree) {
      if (elem.record == 1) {
        elem.record = 2;
        elem.tree->refine();
        for (Tree* child : elem.tree->children()) {
          add_elem(is_deformed, *child).record = 0;
        }
      }
    }
  }
}

void Accessible_mesh::update(std::function<bool(Element&)> refine_criterion, std::function<bool(Element&)> unrefine_criterion)
{
  HEXED_ASSERT(tree, "need a tree to refine");
  auto& elems = elements();
  // decide which elements to (un)refine
  #pragma omp parallel for // parallelize this part since `predicate` could be expensive
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    bool ref = refine_criterion(elem);
    bool unref = unrefine_criterion(elem);
    // record legend: 0 => do nothing; 1 => refine; -1 => unrefine; 2=> already (un)refined (i.e. dead)
    if (ref && !unref) elem.record = 1;
    else if (unref && !ref) elem.record = -1;
    else elem.record = 0;
  }
  int n_orig [2];
  // add new elements
  for (bool is_deformed : {0, 1}) n_orig[is_deformed] = container(is_deformed).element_view().size(); // count how many elements there are before adding, so we know where the new ones start
  for (bool is_deformed : {0, 1}) refine_by_record(is_deformed, 0, n_orig[is_deformed]);
  // ref level smoothing: refine elements with neighbors of ref level more than one level greater than their own
  bool changed;
  do {
    int nd = params.n_dim;
    changed = false;
    for (bool is_deformed : {0, 1}) {
      auto& cont = container(is_deformed);
      auto& cont_elems = cont.element_view();
      for (int i_elem = n_orig[is_deformed]; i_elem < cont_elems.size(); ++i_elem) { // only need to worry about new elements
        auto& elem = cont_elems[i_elem];
        if (elem.record == 0) {
          for (int i_dim = 0; i_dim < nd; ++i_dim) {
            for (int sign = 0; sign < 2; ++sign) {
              Eigen::VectorXi direction = Eigen::VectorXi::Zero(nd);
              direction(i_dim) = math::sign(sign);
              Tree* neighbor = elem.tree->find_neighbor(direction);
              if (neighbor) {
                if (neighbor->refinement_level() < elem.refinement_level() - 1) {
                  if (neighbor->elem) {
                    changed = true;
                    neighbor->elem->record = 1;
                  }
                }
              }
            }
          }
        }
      }
    }
    for (bool is_deformed : {0, 1}) refine_by_record(is_deformed, 0, container(is_deformed).element_view().size());
  } while (changed);
  do {
    int nd = params.n_dim;
    changed = false;
    for (bool is_deformed : {0, 1}) {
      auto& cont = container(is_deformed);
      auto& cont_elems = cont.element_view();
      for (int i_elem = 0; i_elem < n_orig[is_deformed]; ++i_elem) { // FIXME iteration won't quite work right
        auto& elem = cont_elems[i_elem];
        if (elem.record == -1 && elem.tree) {
          bool is_min = true;
          for (int i_dim = 0; i_dim < nd; ++i_dim) is_min = is_min && elem.nominal_position()[i_dim]%2 == 0;
          if (is_min) {
            bool unref = false;
            bool is_def = false; // whether the putative unrefined element will be deformed
            Tree* parent;
            if (!elem.tree->is_root()) {
              unref = true;
              parent = elem.tree->parent();
              for (Tree* child : parent->children()) {
                if (!child->elem) unref = false; // note this will happen if child is not a leaf
                else {
                  unref = unref && child->elem->record == -1;
                  is_def = is_def || child->elem->get_is_deformed();
                }
              }
              if (unref) {
                for (int i_dim = 0; i_dim < nd; ++i_dim) {
                  for (int sign = 0; sign < 2; ++sign) {
                    Eigen::VectorXi direction = Eigen::VectorXi::Zero(nd);
                    direction(i_dim) = math::sign(sign);
                    unref = unref && parent->find_neighbors(direction).size() <= std::size_t(math::pow(2, nd - 1));
                  }
                }
              } else elem.record = 0;
            }
            if (unref) {
              changed = true;
              for (Tree* child : parent->children()) child->elem->record = 2;
              parent->unrefine();
              add_elem(is_def, *parent).record = 0;
            }
          }
        }
      }
    }
  } while (changed);
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    if (elems[i_elem].record == -1) elems[i_elem].record = 0;
  }
  // delete connections to old elements (has to happen before deleting elements or else use after free)
  car.purge_connections();
  def.purge_connections();
  // delete old elements
  n_orig[0] -= car.elems.purge();
  n_orig[1] -= def.elems.purge();
  // connect new elements
  connect_new<         Element>(n_orig[0]);
  connect_new<Deformed_element>(n_orig[1]);
}

}
