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
  matcher_v{car.hanging_vertex_matchers(), def.hanging_vertex_matchers()},
  surf_geom{nullptr}
{
  def.face_con_v = def_face_cons;
}

Accessible_mesh::~Accessible_mesh()
{
  // delete connections before anything else so that deleting elements doesn't create dangling references
  car.purge_connections(always);
  def.purge_connections(always);
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
  car.cons.emplace_back(new Element_face_connection<Element>(el_ar, direction));
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
  def.cons.emplace_back(new Element_face_connection<Deformed_element>(el_ar, direction));
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
  if (is_deformed) {
    def.bound_cons.emplace_back(new Typed_bound_connection<Deformed_element>(def.elems.at(ref_level, element_serial_n), i_dim, face_sign, bc_serial_n));
  } else {
    car.bound_cons.emplace_back(new Typed_bound_connection<Element>(car.elems.at(ref_level, element_serial_n), i_dim, face_sign, bc_serial_n));
  }
}

void Accessible_mesh::disconnect_boundary(int bc_sn)
{
  erase_if(car.bound_cons, [bc_sn](std::unique_ptr<Typed_bound_connection<         Element>>& con){return con->bound_cond_serial_n() == bc_sn;});
  erase_if(def.bound_cons, [bc_sn](std::unique_ptr<Typed_bound_connection<Deformed_element>>& con){return con->bound_cond_serial_n() == bc_sn;});
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
  int n_redundant = 0;
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    for (int i_face = 0; i_face < n_faces; ++i_face) {
      int rec = elems[i_elem].face_record[i_face];
      if (rec == 0) ++n_missing;
      if (rec > 1) n_redundant += rec - 1;
    }
  }
  return {n_redundant, n_missing};
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

void request_connection(Element& elem, int n_dim, int i_dim, bool i_sign, int j_dim, bool j_sign)
{
  // record data at vertex which is on the face to be connected, on the face which was extruded from,
  // and if applicable has the minimum index to satisfy the above consitions.
  int i_vert =   j_sign*math::pow(2, n_dim - 1 - j_dim)
               + (1 - i_sign)*math::pow(2, n_dim - 1 - i_dim);
  auto& record = elem.vertex(i_vert).record;
  // which element it is
  record.push_back(elem.refinement_level());
  record.push_back(elem.record); // `elem.record` = serial number
  // which face needs to be connected
  record.push_back(2*j_dim + j_sign);
  // extrusion direction for deciding which face connections are valid
  record.push_back(2*i_dim + i_sign);
}

void Accessible_mesh::extrude(bool collapse, double offset)
{
  erase_if(vert_ptrs, &Vertex::Non_transferable_ptr::is_null);
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
    // no face has more than one connection (I really hope!) so use numbers greater than one to indentify boundary conditions
    con.element().face_record[2*con.i_dim() + con.inside_face_sign()] = 2 + con.bound_cond_serial_n();
  }

  // request connections with existing extruded elements
  if (tree) {
    def.elems.write_sns();
    for (auto con : extrude_cons) {
      if (con->element(1).tree) {
        auto dir = con->direction();
        for (int j_dim = dir.i_dim[1] + 1; j_dim%nd != dir.i_dim[1]; ++j_dim) {
          j_dim = j_dim%nd;
          for (int face_sign = 0; face_sign < 2; ++face_sign) {
            if (con->element(0).face_record[2*j_dim + face_sign] == 0) {
              request_connection(con->element(0), nd, dir.i_dim[1], dir.face_sign[1], j_dim, face_sign);
            }
          }
        }
      }
    }
  }

  // decide which faces to extrude from
  std::vector<Empty_face> empty_faces;
  auto& elems = def.elements();
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    if (elems[i_elem].tree || !tree) {
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
    elem.record = sn;
    if (collapse) {
      int stride = math::pow(2, nd - 1 - face.i_dim);
      for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
        int i_collapse = i_vert + (face.face_sign - (i_vert/stride)%2)*stride;
        elem.vertex(i_vert).pos = face.elem.vertex(i_collapse).pos;
      }
    }
    std::array<Deformed_element*, 2> el_arr {&elem, &face.elem};
    def.cons.emplace_back(new Element_face_connection<Deformed_element>(el_arr, dir));
    extrude_cons.push_back(def.cons.back().get());
    // record the faces that still need to be connected at a vertex which is guaranteed to be shared with prospective neighbors
    for (int j_dim = face.i_dim + 1; j_dim%nd != face.i_dim; ++j_dim) {
      j_dim = j_dim%nd;
      for (int face_sign = 0; face_sign < 2; ++face_sign) {
        int face_rec = face.elem.face_record[2*j_dim + face_sign];
        if (face_rec >= 2) {
          // if parent element has boundary connections on other faces
          def.bound_cons.emplace_back(new Typed_bound_connection<Deformed_element>(elem, j_dim, face_sign, face_rec - 2));
        } else request_connection(elem, nd, face.i_dim, face.face_sign, j_dim, face_sign);
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
    #if 0
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
      HEXED_ASSERT(verts[i_vert].record.empty(), format_str(100, "vertex detected with %lu unprocessed connection requests",
                                                            verts[i_vert].record.size()/n_record).c_str());
    }
    #endif
  }
  for (auto con_plan : con_plans) {
    connect_deformed(con_plan.ref_level, con_plan.serial_ns, con_plan.dir);
  }
}

void Accessible_mesh::connect_rest(int bc_sn)
{
  auto& elem_seq = elements();
  // locate unconnected faces
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elem_seq.size(); ++i_elem) {
    for (int i_face = 0; i_face < 2*params.n_dim; ++i_face) {
      elem_seq[i_elem].face_record[i_face] = 0;
    }
  }
  car.record_connections();
  def.record_connections();
  // make connections
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
  if (is_deformed) {
    t.def_elem = &def.elems.at(t.refinement_level(), sn);
  }
  return elem;
}

void Accessible_mesh::add_tree(std::vector<Flow_bc*> extremal_bcs)
{
  // take ownership of bcs (do this first to avoid memory leak)
  std::vector<int> new_tree_bcs;
  //! \todo this could, in theory, be a resource leak because these are never erased if an exception is thrown...
  for (Flow_bc* fbc : extremal_bcs) new_tree_bcs.push_back(add_boundary_condition(fbc, new Nominal_pos));
  // make assertions
  HEXED_ASSERT(int(extremal_bcs.size()) == 2*params.n_dim, "`extremal_bcs` has wrong number of elements");
  HEXED_ASSERT(!tree, "each `Mesh` may only contain one tree");
  // add the tree
  tree_bcs = new_tree_bcs;
  tree.reset(new Tree(params.n_dim, root_sz));
  auto& elem = add_elem(false, *tree);
  int sn = elem.record;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      connect_boundary(0, 0, sn, i_dim, sign, tree_bcs[2*i_dim + sign]);
    }
  }
}

bool Accessible_mesh::is_surface(Tree* t)
{
  if (!surf_geom) return false;
  Mat<> center = t->nominal_position() + t->nominal_size()/2*Mat<>::Ones(params.n_dim);
  double dist = (center - surf_geom->nearest_point(center)).norm();
  return dist < std::sqrt(params.n_dim)/2*t->nominal_size();
}

void Accessible_mesh::set_surface(Surface_geom* geometry, Flow_bc* surface_bc, Eigen::VectorXd flood_fill_start)
{
  // take ownership of the surface geometries (do this first to avoid memory leak)
  surf_bc_sn = add_boundary_condition(surface_bc, new Geom_mbc(geometry));
  surf_geom = geometry;
  if (!tree) return;
  // identify surface elements
  auto& elems = elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    if (elem.tree) if (is_surface(elem.tree)) elem.tree->set_status(0);
  }
  // pefrorm flood fill
  Tree* start = tree->find_leaf(flood_fill_start);
  if (!start) start = tree.get();
  start->flood_fill(1);
  // delete stuff
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    elem.record = 0;
    if (elem.tree) if (elem.tree->get_status() != 1) elem.record = 2;
  }
  deform();
  purge();
  connect_new<         Element>(0);
  connect_new<Deformed_element>(0);
  extrude();
  connect_rest(surf_bc_sn);
}

// forms connections for new tree elements of type `element_t` starting with the `starting_at`th one
template<typename element_t>
void Accessible_mesh::connect_new(int start_at)
{
  auto& m = mbt<element_t>();
  auto elems = m.elems.elements();
  int nd = params.n_dim;
  // helper function for connecting refined elements
  auto connect_refined = [&](Element& elem, int i_dim, int sign, std::vector<Tree*> neighbors) {
    HEXED_ASSERT(int(neighbors.size()) == math::pow(2, nd - 1), format_str(100, "bad number of neighbors %lu (thanks for nothing, ref level smoother)", neighbors.size()))
    bool is_def = elem.get_is_deformed();
    for (Tree* neighbor : neighbors) {
      HEXED_ASSERT(neighbor->elem, "hanging-node connection with nonexistant elements");
      is_def = is_def && neighbor->elem->get_is_deformed();
    }
    if (is_def) {
      std::vector<Deformed_element*> fine;
      for (Tree* neighbor : neighbors) fine.push_back(neighbor->def_elem);
      def.ref_face_cons[nd - 1].emplace_back(new Refined_connection<Deformed_element>(elem.tree->def_elem, fine, Con_dir<Deformed_element>{{i_dim, i_dim}, {!sign, bool(sign)}}));
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
            // if this element is at the boundary of the tree (as opposed to a surface geometry boundary) set an extremal boundary condition
            if (neighbors.empty()) m.bound_cons.emplace_back(new Typed_bound_connection<element_t>(elem, i_dim, sign, tree_bcs[2*i_dim + sign]));
            // otherwise, if the element has not only a tree neighbor but also an element neighbor...
            else if (neighbors[0]->elem) {
              if (neighbors.size() == 1) {
                auto& other = *neighbors[0]->elem;
                // if ref levels are the same, make a conformal connection
                if (other.refinement_level() == elem.refinement_level()) {
                  if (elem.get_is_deformed() && other.get_is_deformed()) {
                    std::array<Deformed_element*, 2> el_ar {elem.tree->def_elem, neighbors[0]->def_elem};
                    def.cons.emplace_back(new Element_face_connection<Deformed_element>(el_ar, Con_dir<Deformed_element>{{i_dim, i_dim}, {bool(sign), !sign}}));
                  } else {
                    std::array<Element*, 2> el_ar;
                    el_ar[!sign] = &elem;
                    el_ar[sign] = &other;
                    car.cons.emplace_back(new Element_face_connection<Element>(el_ar, Con_dir<Element>{i_dim}));
                  }
                } else {
                  // if neighbor is coarser, form a hanging node connection
                  // but only if this is the fine element with the lowest coordinates, to prevent redundant connections from all the fine elements
                  bool is_min_corner = true;
                  for (int j_dim = 0; j_dim < nd; ++j_dim) if (j_dim != i_dim) {
                    is_min_corner = is_min_corner && elem.tree->coordinates()[j_dim]%2 == 0;
                  }
                  if (is_min_corner) connect_refined(other, i_dim, sign, other.tree->find_neighbors(-direction));
                }
              } else {
                // if neighbors are finer, form a hanging node connection
                connect_refined(elem, i_dim, !sign, neighbors);
              }
            }
          }
        }
      }
    }
  }
}

// performs the actual refinement for all elements where the record has been set to 1
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

/* returns:
 * - do any neighbors exist?
 * - do all neighbors exist?
 * - are any of the neighbors more than 1 ref level higher than `t`?
 */
std::array<int, 3> eval_neigbors(Tree* t, Eigen::VectorXi direction)
{
  std::array<int, 3> result {false, false, false};
  auto neighbors = t->find_neighbors(direction);
  if (!neighbors.empty()) {
    result[1] = true;
    for (Tree* n : neighbors) {
      bool exists = false;
      if (n->elem) if (n->elem->record != 2) {
        if (n->refinement_level() > t->refinement_level() + 1) result[2] = true; // if ref level difference is greater than 1, need to refine
        exists = true;
      }
      result[0] = result[0] || exists;
      result[1] = result[1] && exists;
    }
  }
  return result;
}

// does this tree element need to be refined to satisfy ref level smoothness
bool Accessible_mesh::needs_refine(Tree* t)
{
  for (int i_face = 0; i_face < 2*params.n_dim; ++i_face) {
    auto eval = eval_neigbors(t, math::direction(params.n_dim, i_face));
    if (eval[2]) return true; // if there is a face neighbor with ref level more than 1 greater, need to refine
    if (t->elem) {
      if (eval[0] && !eval[1]) return true; // if there is a partially-exposed face, need to refine
      // if this face is exposed, also check the diagonal neighbors since they could be extrusion neighbors
      else if (!eval[1]) {
        for (int j_dim = 0; j_dim < params.n_dim; ++j_dim) if (j_dim != i_face/2) {
          for (int face_sign = 0; face_sign < 2; ++face_sign) {
            auto dir = math::direction(params.n_dim, i_face);
            dir(j_dim) = math::sign(face_sign);
            auto eval1 = eval_neigbors(t, dir);
            // again, ref level difference > 1 or partially exposed -> refine
            if (eval1[2] || (eval1[0] && !eval1[1])) return true;
          }
        }
      }
    }
  }
  // otherwise, we're good
  return false;
}

bool has_existent_children(Tree* t)
{
  bool has = t->elem;
  for (Tree* child : t->children()) has = has || has_existent_children(child);
  return has;
}

// whether an element is currently set to be deformed at the end of the `update` sweep
bool is_def(Element& elem)
{
  return (elem.get_is_deformed() && elem.record == 0) || (!elem.get_is_deformed() && elem.record == 3);
}

void Accessible_mesh::deform()
{
  int nd = params.n_dim;
  auto& elems = elements();
  // start with all elements as cartesian
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    if (elem.record != 2 && elem.get_is_deformed() && elem.tree) elem.record = 3;
  }
  // deform all boundary elements and certain of their face neighbors
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    if (elem.record != 2 && elem.tree) {
      // figure out if this element has a face on the boundary
      std::vector<Tree*> neighbors [6];
      bool surface [6] {};
      bool boundary = false;
      for (int i_face = 0; i_face < 2*nd; ++i_face) {
        neighbors[i_face] = elem.tree->find_neighbors(math::direction(nd, i_face));
        if (!neighbors[i_face].empty()) {
          if (neighbors[i_face][0]->elem) surface[i_face] = neighbors[i_face][0]->elem->record == 2;
          else surface[i_face] = true;
        }
        boundary = boundary || surface[i_face];
      }
      // if it does, deform it and all face neighbors which are not opposite the boundary face
      if (boundary) {
        #pragma omp atomic write
        elem.record = 3*!elem.get_is_deformed();
        for (int i_face = 0; i_face < 2*nd; ++i_face) {
          bool def_face = false;
          for (int j_face = 0; j_face < 2*nd; ++j_face) {
            if (j_face/2 != i_face/2 && surface[j_face]) def_face = true;
          }
          if (def_face) {
            for (Tree* n : neighbors[i_face]) if (n->elem) if (n->elem->record != 2) {
              #pragma omp atomic write
              n->elem->record = 3*!n->elem->get_is_deformed();
            }
          }
        }
      }
    }
  }
  // if any refined faces have some elements cartesian and some deformed, make them all deformed
  bool changed;
  do {
    changed = false;
    #pragma omp parallel for reduction(||:changed)
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      if (elem.tree && is_def(elem)) {
        for (int i_face = 0; i_face < 2*nd; ++i_face) {
          auto neighbors = elem.tree->find_neighbors(math::direction(nd, i_face));
          bool any_deformed = false;
          bool all_deformed = true;
          for (Tree* neighbor : neighbors) {
            if (neighbor->elem) {
              any_deformed = any_deformed || is_def(*neighbor->elem);
              all_deformed = all_deformed && is_def(*neighbor->elem);
            }
          }
          if (any_deformed && !all_deformed) {
            for (Tree* neighbor : neighbors) {
              if (neighbor->elem) {
                changed = true;
                #pragma omp atomic write
                neighbor->elem->record = 3*!neighbor->elem->get_is_deformed();
              }
            }
          }
        }
      }
    }
  } while (changed);
  // add new elements
  for (bool is_deformed : {0, 1}) {
    auto& cont = container(is_deformed);
    auto& cont_elems = cont.element_view();
    int sz = cont_elems.size();
    for (int i_elem = 0; i_elem < sz; ++i_elem) {
      auto& elem = cont_elems[i_elem];
      if (elem.record == 3) {
        add_elem(!is_deformed, *elem.tree).record = 0;
        elem.record = 2;
        elem.tree = nullptr;
      }
    }
  }
}

void Accessible_mesh::purge()
{
  // delete obsolete elements of `extrude_cons`
  erase_if(extrude_cons, [](Element_face_connection<Deformed_element>* con){return con->element(0).record == 2 || con->element(1).record == 2;});
  // delete connections to old elements (has to happen before deleting elements or else use after free)
  car.purge_connections();
  def.purge_connections();
  // delete old elements
  car.elems.purge();
  def.elems.purge();
  // delete dangling vertex pointers
  erase_if(vert_ptrs, &Vertex::Non_transferable_ptr::is_null);
}

void Accessible_mesh::update(std::function<bool(Element&)> refine_criterion, std::function<bool(Element&)> unrefine_criterion)
{
  /* `Element::record` is used to identify which elements are going to be modified.
   * 0 => do nothing
   * 1 => refine
   * -1 => unrefine
   * 2 => delete
   * 3 => toggle deformity
   */
  HEXED_ASSERT(tree, "need a tree to refine");
  int nd = params.n_dim;
  auto& elems = elements();
  // decide which elements to (un)refine
  #pragma omp parallel for // parallelize this part since `predicate` could be expensive
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    elem.record = 0;
    if (elem.tree) {
      bool ref = refine_criterion(elem);
      bool unref = unrefine_criterion(elem);
      if (ref && !unref) elem.record = 1;
      else if (unref && !ref) elem.record = -1;
    }
  }
  int n_orig [2];
  // refine elements
  for (bool is_deformed : {0, 1}) n_orig[is_deformed] = container(is_deformed).element_view().size(); // count how many elements there are before adding, so we know where the new ones start
  for (bool is_deformed : {0, 1}) refine_by_record(is_deformed, 0, n_orig[is_deformed]);
  // synchronize refinement level of surface elements with their non-surface neighbors in preparation for incremental flood fill
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) { // only need to worry about new elements
    auto& elem = elems[i_elem];
    if (elem.tree) {
      for (int i_face = 0; i_face < 2*nd; ++i_face) {
        Tree* neighbor = elem.tree->find_neighbor(math::direction(nd, i_face));
        if (neighbor) {
          if (!neighbor->elem) {
            if (neighbor->refinement_level() > elem.refinement_level()) {
              Tree* p = neighbor->parent();
              bool can_unref = true;
              for (Tree* child : p->children()) can_unref = can_unref && !child->elem;
              if (can_unref && !needs_refine(p)) p->unrefine();
            } else if (neighbor->refinement_level() < elem.refinement_level() - 1) neighbor->refine();
            else if (neighbor->refinement_level() < elem.refinement_level()) {
              int min_rl = std::numeric_limits<int>::max();
              for (int j_face = 0; j_face < 2*nd; ++j_face) {
                Tree* n = neighbor->find_neighbor(math::direction(nd, j_face));
                if (n) if (n->elem) min_rl = std::min(min_rl, n->refinement_level());
              }
              if (neighbor->refinement_level() < min_rl) neighbor->refine();
            }
          }
        }
      }
    }
  }
  bool changed;
  // incremental flood fill
  do {
    changed = false;
    for (bool is_deformed : {0, 1}) {
      auto& cont = container(is_deformed);
      auto& cont_elems = cont.element_view();
      int sz = cont_elems.size();
      for (int i_elem = 0; i_elem < sz; ++i_elem) {
        auto& elem = cont_elems[i_elem];
        if (elem.tree && elem.record == 0) {
          for (int i_face = 0; i_face < 2*nd; ++i_face) {
            for (Tree* neighbor : elem.tree->find_neighbors(math::direction(nd, i_face))) {
              if (!neighbor->elem) if (!is_surface(neighbor)) {
                changed = true;
                add_elem(is_deformed, *neighbor).record = 0;
              }
            }
          }
        }
      }
    }
  } while (changed);
  // unrefinement
  do {
    changed = false;
    for (bool is_deformed : {0, 1}) {
      auto& cont_elems = container(is_deformed).element_view();
      for (int i_elem = 0; i_elem < n_orig[is_deformed]; ++i_elem) {
        auto& elem = cont_elems[i_elem];
        if (elem.record == -1 && elem.tree) {
          bool unref = false;
          bool is_def = false; // whether the putative unrefined element will be deformed
          Tree* parent;
          if (!elem.tree->is_root()) { // can't unrefine the root
            unref = true;
            parent = elem.tree->parent();
            // only unrefine if all the existing siblings agree and have the same ref level
            for (Tree* child : parent->children()) {
              if (!child->is_leaf() && has_existent_children(child)) unref = false;
              else if (child->elem) if (child->elem->record != 2) {
                unref = unref && child->elem->record == -1;
                is_def = is_def || child->elem->get_is_deformed();
              }
            }
            if (unref) unref = unref && !needs_refine(parent); // don't unrefine if it would violate ref level smoothness
            else elem.record = 0; // if we didn't unrefine because of one of the siblings, set the record to 0 to avoid redundant checks in future sweeps
          }
          // perform unrefinement
          if (unref) {
            changed = true;
            for (Tree* child : parent->children()) {
              if (child->elem) {
                child->elem->record = 2;
                child->elem->tree = nullptr;
              }
            }
            parent->unrefine();
            add_elem(is_def, *parent).record = 0;
          }
        }
      }
    }
  } while (changed);
  // set the record straight for any elements denied refinement
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    if (elems[i_elem].record == -1) elems[i_elem].record = 0;
  }
  // un-flood-fill any new surface elements
  for (bool is_deformed : {0, 1}) {
    auto& cont_elems = container(is_deformed).element_view();
    #pragma omp parallel for
    for (int i_elem = n_orig[is_deformed]; i_elem < cont_elems.size(); ++i_elem) {
      auto& elem = cont_elems[i_elem];
      if (elem.tree && elem.record != 2) if (is_surface(elem.tree)) elem.record = 2;
    }
  }
  // ref level smoothing: refine elements to satisfy solver requirements on neighbors
  do {
    changed = false;
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      if (elem.record != 2 && elem.tree) {
        if (needs_refine(elem.tree)) {
          changed = true;
          elem.record = 1;
        }
      }
    }
    for (bool is_deformed : {0, 1}) refine_by_record(is_deformed, 0, container(is_deformed).element_view().size());
  } while (changed);
  deform();
  // delete extruded elements
  #pragma omp parallel for
  for (auto con : extrude_cons) {
    bool del = false;
    if (con->element(1).record == 2) del = true;
    else if (con->element(1).tree) {
      Tree* neighbor = con->element(1).tree->find_neighbor(math::direction(nd, con->direction().i_face(1)));
      if (neighbor) if (neighbor->elem) if (neighbor->elem->record != 2) del = true;
    }
    if (del) con->element(0).record = 2;
  }
  for (bool is_deformed : {0, 1}) {
    auto& cont = container(is_deformed);
    auto& cont_elems = cont.element_view();
    int deleted = 0;
    #pragma omp parallel for reduction(+:deleted)
    for (int i_elem = 0; i_elem < n_orig[is_deformed]; ++i_elem) {
      deleted += cont_elems[i_elem].record == 2;
    }
    n_orig[is_deformed] -= deleted;
  }
  purge();
  // connect new elements
  connect_new<         Element>(n_orig[0]);
  connect_new<Deformed_element>(n_orig[1]);
  extrude();
  connect_rest(surf_bc_sn);
}

}
