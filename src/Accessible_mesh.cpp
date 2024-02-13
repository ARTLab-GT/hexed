#include <Accessible_mesh.hpp>
#include <math.hpp>
#include <Row_index.hpp>
#include <erase_if.hpp>
#include <utils.hpp>
#include <Gauss_legendre.hpp>
#include <H5Cpp.h>
#include <filesystem>
#include <fstream>

namespace hexed
{

Element_container& Accessible_mesh::container(bool is_deformed)
{
  Element_container* containers [] {&car.elems, &def.elems};
  return *containers[is_deformed];
}

template<> Mesh_by_type<         Element>& Accessible_mesh::mbt() {return car;}
template<> Mesh_by_type<Deformed_element>& Accessible_mesh::mbt() {return def;}

void Accessible_mesh::id_boundary_verts()
{
  auto verts = vertices();
  #pragma omp parallel for
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
    verts[i_vert].record.clear();
  }
  for (auto& b_verts : boundary_verts) {
    erase_if(b_verts, [](Vertex::Non_transferable_ptr& ptr) {
      if (!ptr) return true;
      return !ptr->needs_smooth();});
  }
  for (unsigned bc_sn = 0; bc_sn < boundary_verts.size(); ++bc_sn) {
    #pragma omp parallel for
    for (auto& vert : boundary_verts[bc_sn]) {
      vert->record.push_back(bc_sn);
    }
    auto& cons = def.boundary_connections();
    for (int i_con = 0; i_con < cons.size(); ++i_con) {
      auto& con = cons[i_con];
      if (con.bound_cond_serial_n() == int(bc_sn)) {
        for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
          if ((i_vert/math::pow(2, params.n_dim - 1 - con.i_dim()))%2 == con.inside_face_sign()) {
            auto& vert = con.element().vertex(i_vert);
            Lock::Acquire a(vert.lock);
            if (!std::count(vert.record.begin(), vert.record.end(), bc_sn)) {
              vert.record.push_back(bc_sn);
              if (vert.needs_smooth()) boundary_verts[bc_sn].emplace_back(vert);
            }
          }
        }
      }
    }
  }
}

void Accessible_mesh::snap_vertices()
{
  // this is a terrible way to find the max ref level. future self please fix
  double min_sz = huge;
  #pragma omp parallel for reduction(min : min_sz)
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    min_sz = std::max(min_sz, elems[i_elem].nominal_size());
  }
  if (tree) {
    auto snap_extremes = [this]() {
      for (int i_bc = 0; i_bc < 2*params.n_dim; ++i_bc) {
        int bc_sn = tree_bcs[i_bc];
        #pragma omp parallel for
        for (auto& vert : boundary_verts[bc_sn]) {
          int i_dim = i_bc/2;
          int sign = i_bc%2;
          vert->pos(i_dim) = tree->origin()(i_dim) + sign*tree->nominal_size();
        }
      }
    };
    // snap extremes before and after snapping to surface to ensure exact extreme snapping and accurate surface snapping
    snap_extremes();
    if (surf_geom) {
      // if i just implement a better way to find the distance guess, we won't need all these atomics
      #pragma omp parallel for
      for (auto& vert : boundary_verts[surf_bc_sn]) {
        Mat<> pos(params.n_dim);
        for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
          #pragma omp atomic read
          pos(i_dim) = vert->pos(i_dim);
        }
        double dist_guess = min_sz;
        for (const Vertex& neighb : vert->get_neighbors()) {
          Mat<> neighb_pos(params.n_dim);
          for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
            #pragma omp atomic read
            neighb_pos(i_dim) = neighb.pos(i_dim);
          }
          dist_guess = std::max(dist_guess, (neighb_pos - pos).norm());
        }
        pos = surf_geom->nearest_point(pos, huge, dist_guess).point();
        for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
          #pragma omp atomic write
          vert->pos(i_dim) = pos(i_dim);
        }
      }
    }
    snap_extremes();
  } else {
    auto& bc_cons {boundary_connections()};
    #pragma omp parallel for
    for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
      int bc_sn = bc_cons[i_con].bound_cond_serial_n();
      boundary_condition(bc_sn).mesh_bc->snap_vertices(bc_cons[i_con]);
    }
  }
  // vertex relaxation/snapping will cause hanging vertices to drift away from hanging vertex faces they are supposed to be coincident with
  // so now we put them back where they belong
  auto& matchers = hanging_vertex_matchers();
  #pragma omp parallel for
  for (int i_match = 0; i_match < matchers.size(); ++i_match) {
    matchers[i_match].match(&Element::vertex_position<0>);
    matchers[i_match].match(&Element::vertex_position<1>);
    matchers[i_match].match(&Element::vertex_position<2>);
  }
}

Accessible_mesh::Accessible_mesh(Storage_params params_arg, double root_size_arg) :
  params{params_arg},
  n_vert{math::pow(2, params.n_dim)},
  root_sz{root_size_arg},
  car{params, root_sz},
  def{params, root_sz},
  def_as_car{def.elements()},
  elems{car.elements(), def_as_car},
  kernel_elems{elems},
  elem_cons{car.element_connections(), def.element_connections()},
  bound_face_cons{car.bound_face_con_view, def.bound_face_con_view},
  bound_cons{car.boundary_connections(), def.boundary_connections()},
  def_face_cons{def.elem_face_con_v, bound_face_cons},
  ref_face_v{car.refined_faces(), def.refined_faces()},
  matcher_v{car.hanging_vertex_matchers(), def.hanging_vertex_matchers()},
  surf_bc_sn{-1}, // set to -1 to prevent uninitialized comparisons
  surf_geom{nullptr},
  verts_are_reset{false},
  buffer_dist{std::sqrt(params.n_dim)/2}
{
  def.face_con_v = def_face_cons;
}

Accessible_mesh::~Accessible_mesh()
{
  // delete connections before anything else so that deleting elements doesn't create dangling references
  car.purge_connections(criteria::always);
  def.purge_connections(criteria::always);
}

int Accessible_mesh::add_element(int ref_level, bool is_deformed, std::vector<int> position, Mat<> origin)
{
  int sn = container(is_deformed).emplace(ref_level, position, origin);
  Element& elem = element(ref_level, is_deformed, sn);
  for (int i_vert = 0; i_vert < n_vert; ++i_vert) vert_ptrs.emplace_back(elem.vertex(i_vert));
  return sn;
}

int Accessible_mesh::add_element(int ref_level, bool is_deformed, std::vector<int> position)
{
  return add_element(ref_level, is_deformed, position, Mat<>::Zero(params.n_dim));
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
  boundary_verts.emplace_back();
  bound_conds.push_back({std::unique_ptr<Flow_bc>{flow_bc}, std::unique_ptr<Mesh_bc>{mesh_bc}});
  // no reason to delete boundary conditions, so the serial number can just be the index
  return bound_conds.size() - 1;
}

void Accessible_mesh::connect_boundary(int ref_level, bool is_deformed, int element_serial_n, int i_dim, int face_sign, int bc_serial_n)
{
  // create boundary condition
  HEXED_ASSERT(bc_serial_n < int(bound_conds.size()), "demand for non-existent `Boundary_condition`");
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

void Accessible_mesh::cleanup()
{
  purge();
  id_smooth_verts();
  id_boundary_verts();
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
  HEXED_ASSERT(!verts_are_reset, "looks like you have a `Mesh::Reset_vertices` you forgot to destroy");
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

void Accessible_mesh::extrude(bool collapse, double offset, bool force)
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
  {
    auto verts = vertices();
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
      HEXED_ASSERT(verts[i_vert].record.size()%4 == 0, "vertex has wrong number of records");
    }
  }

  // decide which faces to extrude from
  std::vector<Empty_face> empty_faces;
  auto& elems = def.elements();
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    if (elems[i_elem].tree || (!tree || force)) {
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
    int sn = add_element(ref_level, true, nom_pos, face.elem.origin);
    Con_dir<Deformed_element> dir {{face.i_dim, face.i_dim}, {!face.face_sign, bool(face.face_sign)}};
    auto& elem = def.elems.at(ref_level, sn);
    elem.record = sn;
    elem.needs_snapping = !force;
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
    if (offset > 0) {
      // readjust face node adjustments to account for offset
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
      double* state [] {face.elem.state(), elem.state()};
      // interpolate data from original element to new ones
      for (int i_elem : {1, 0}) { // iterate in reverse order since new states for both elements depend on element 0
        Gauss_legendre basis(params.row_size); //! \todo apparently the mesh needs to know about the basis after all...
        double width = 1 - i_elem + math::sign(i_elem)*offset;
        Mat<dyn, dyn> interp = basis.interpolate(basis.nodes()*width + Mat<>::Constant(params.row_size, (i_elem == face.face_sign)*(1 - width)));
        for (Row_index index(nd, params.row_size, face.i_dim); index; ++index) {
          Eigen::Map<Mat<dyn, dyn>, 0, Eigen::Stride<dyn, dyn>> row_read (state[0     ] + index.i_qpoint(0), params.row_size, params.n_var_numeric(), Eigen::Stride<dyn, dyn>(params.n_qpoint(), index.stride));
          Eigen::Map<Mat<dyn, dyn>, 0, Eigen::Stride<dyn, dyn>> row_write(state[i_elem] + index.i_qpoint(0), params.row_size, params.n_var_numeric(), Eigen::Stride<dyn, dyn>(params.n_qpoint(), index.stride));
          row_write = interp*row_read;
        }
      }
    }
  }
  {
    auto verts = vertices();
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
      HEXED_ASSERT(verts[i_vert].record.size()%4 == 0, "vertex has wrong number of records");
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
  {
    auto verts = vertices();
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
      HEXED_ASSERT(verts[i_vert].record.size()%4 == 0, "vertex has wrong number of records");
    }
  }

  // plan connections to make (don't make them yet, because that could result in `eat`ing vertices which have not been visited,
  // and ultimately dereferencing null pointers)
  std::vector<Connection_plan> con_plans;
  std::vector<Refined_connection_plan> ref_con_plans;

  #define VERTEX_LOOP(code) \
    /* connections without hanging nodes */ \
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert) { \
      auto& vert {verts[i_vert]}; \
      /* first make connections where dimension matches and then make connections among differing dimensions. */ \
      /* This order prevents incorrect connections where both same-dimension and different-dimension candidates are available. */ \
      for (bool (*aligned)(Con_dir<Deformed_element>, std::array<int, 2>) : {&aligned_same_dim, &aligned_different_dim}) { \
        /* iterate through every possible pair of records created by an extruded elements above */ \
        for (int i_record = 0; i_record < int(vert.record.size()); i_record += n_record) { \
          for (int j_record = i_record + n_record; j_record < int(vert.record.size()); j_record += n_record) { \
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
        HEXED_ASSERT(fine_vert.record.size() == n_record, \
                     format_str(1000, "`fine_vert.record.size() == %li != n_record == %li` (position = [%e, %e, %e])", \
                                fine_vert.record.size(), n_record, \
                                fine_vert.pos[0], fine_vert.pos[1], fine_vert.pos[2]).c_str()); \
        sn[1] = fine_vert.record[1]; \
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
  int sn = add_element(t.refinement_level(), is_deformed, std::vector<int>(np.begin(), np.end()), t.origin());
  auto& elem = element(t.refinement_level(), is_deformed, sn);
  elem.record = sn; // put the serial number in the record so it can be used for connections
  elem.tree.pair(t.elem);
  if (is_deformed) {
    t.def_elem = &def.elems.at(t.refinement_level(), sn);
  }
  return elem;
}

void Accessible_mesh::create_tree(std::vector<Flow_bc*> extremal_bcs, Mat<> origin)
{
  // take ownership of bcs (do this first to avoid memory leak)
  std::vector<int> new_tree_bcs;
  //! \todo this could, in theory, be a resource leak because these are never erased if an exception is thrown...
  for (Flow_bc* fbc : extremal_bcs) new_tree_bcs.push_back(add_boundary_condition(fbc, new Nominal_pos));
  HEXED_ASSERT(int(extremal_bcs.size()) == 2*params.n_dim, "`extremal_bcs` has wrong number of elements");
  HEXED_ASSERT(!tree, "each `Mesh` may only contain one tree");
  // add the tree
  tree_bcs = new_tree_bcs;
  tree.reset(new Tree(params.n_dim, root_sz, origin));
}

void Accessible_mesh::add_tree(std::vector<Flow_bc*> extremal_bcs, Mat<> origin)
{
  create_tree(extremal_bcs, origin);
  auto& elem = add_elem(false, *tree);
  int sn = elem.record;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      connect_boundary(0, 0, sn, i_dim, sign, tree_bcs[2*i_dim + sign]);
    }
  }
}

bool Accessible_mesh::intersects_surface(Tree* t)
{
  if (!surf_geom) return false;
  Mat<> center = t->nominal_position() + t->nominal_size()/2*Mat<>::Ones(params.n_dim);
  return !surf_geom->nearest_point(center, buffer_dist*t->nominal_size()).empty();
}

bool Accessible_mesh::is_surface(Tree* t)
{
  return t->get_status() == 0;
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
    if (elem.tree) if (intersects_surface(elem.tree)) elem.tree->set_status(0);
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
  delete_bad_extrusions();
  deform();
  purge();
  connect_new<         Element>(0);
  connect_new<Deformed_element>(0);
  extrude(true);
  connect_rest(surf_bc_sn);
  id_boundary_verts();
  snap_vertices();
  id_smooth_verts();
}

void Accessible_mesh::set_unref_locks(std::function<bool(Element&)> lock_if)
{
  auto& elems = elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].unrefinement_locked = lock_if(elems[i_elem]);
  }
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
          if (!elem.is_connected(2*i_dim + sign)) { // if this face hasn't been connected already
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

void Accessible_mesh::refine_set_status(Tree* t)
{
  t->refine();
  for (Tree* child : t->children()) {
    child->set_status(intersects_surface(child) - 1);
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
        refine_set_status(elem.tree);
        for (Tree* child : elem.tree->children()) {
          add_elem(is_deformed, *child).record = 0;
        }
      }
    }
  }
}

bool exists(Tree* tree)
{
  if (tree) if (tree->elem) {
    int record;
    #pragma omp atomic read
    record = tree->elem->record;
    if (record != 2) return true;
  }
  return false;
}

// does this tree element need to be refined to satisfy ref level smoothness
bool Accessible_mesh::needs_refine(Tree* t)
{
  for (int i_face = 0; i_face < 2*params.n_dim; ++i_face) {
    // if there is a face neighbor with ref level more than 1 greater, need to refine
    auto dir = math::direction(params.n_dim, i_face);
    auto neighbors = t->find_neighbors(dir);
    int rl = t->refinement_level();
    auto too_fine = [rl](Tree* ptr){return exists(ptr) && ptr->refinement_level() > rl + 1;};
    if (std::any_of(neighbors.begin(), neighbors.end(), too_fine)) return true;
    if (t->elem) {
      // also check the diagonal neighbors since they could be extrusion neighbors
      for (int j_face = 0; j_face < 2*(i_face/2); ++j_face) {
        auto d = dir;
        d(j_face/2) = math::sign(j_face%2);
        neighbors = t->find_neighbors(d);
        // again, ref level difference > 1 or partially exposed -> refine
        if (std::any_of(neighbors.begin(), neighbors.end(), too_fine)) return true;
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

// delete elements that would create pathological extrusion topology
void Accessible_mesh::delete_bad_extrusions()
{
  int nd = params.n_dim;
  auto& elems = elements();
  bool changed;
  do {
    changed = false;
    #pragma omp parallel for reduction(||:changed)
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      int record;
      #pragma omp atomic read
      record = elem.record;
      if (elem.tree && record != 2) {
        bool exposed [6];
        for (int i_face = 0; i_face < 2*nd; ++i_face) {
          // if any hanging-node face is partially covered by fine elements, delete the fine elements
          auto neighbors = elem.tree->find_neighbors(math::direction(nd, i_face));
          bool all_exist = std::all_of(neighbors.begin(), neighbors.end(), exists);
          exposed[i_face] = !all_exist;
          if (neighbors.size() > 1) {
            if (std::any_of(neighbors.begin(), neighbors.end(), exists) && !all_exist) {
              changed = true;
              for (Tree* n : neighbors) if (n->elem) {
                #pragma omp atomic write
                n->elem->record = 2;
              }
            }
          }
          // for all faces that share an edge with `i_face` (but only the ones with lower index to avoid redundancy)
          for (int j_face = 0; j_face < 2*(i_face/2); ++j_face) {
            if (exposed[i_face] || exposed[j_face]) {
              auto dir = math::direction(nd, i_face);
              dir(j_face/2) = math::sign(j_face%2);
              auto diag_neighbs = elem.tree->find_neighbors(dir);
              // if there are elements that are connected diagonally but have no mutual face neighbors, delete at least one of them
              if (exposed[i_face] && exposed[j_face]) {
                for (Tree* n : diag_neighbs) if (exists(n)) {
                  // To make results repeatable, if the elements have different refinement levels, delete the finer one(s).
                  // If they have the same refinement level, delete both
                  if (n->refinement_level() >= elem.refinement_level()) {
                    changed = true;
                    #pragma omp atomic write
                    n->elem->record = 2;
                    if (n->refinement_level() == elem.refinement_level()) {
                      #pragma omp atomic write
                      elem.record = 2;
                    }
                  }
                }
                for (int k_face = 0; k_face < 2*(j_face/2); ++k_face) if (exposed[k_face]) {
                  auto d = dir;
                  d(k_face/2) = math::sign(k_face%2);
                  auto neighbs = elem.tree->find_neighbors(d);
                  for (Tree* n : neighbs) if (exists(n)) {
                    if (n->refinement_level() >= elem.refinement_level()) {
                      changed = true;
                      #pragma omp atomic write
                      n->elem->record = 2;
                      if (n->refinement_level() == elem.refinement_level()) {
                        #pragma omp atomic write
                        elem.record = 2;
                      }
                    }
                  }
                }
              } else if (nd == 3) {
                // if the edge is partially covered with fine elements, delete them
                if (  !std::all_of(diag_neighbs.begin(), diag_neighbs.end(), exists)
                    && std::any_of(diag_neighbs.begin(), diag_neighbs.end(), exists)) {
                  changed = true;
                  for (Tree* n : diag_neighbs) {
                    if (n->elem) {
                      #pragma omp atomic write
                      n->elem->record = 2;
                    }
                  }
                }
              }
            }
          }
        }
        // delete elements with exposed faces, edges, or vertices that face away from the surface geometry
        if (surf_geom) {
          bool exp = false;
          for (int i_face = 0; i_face < 2*nd; ++i_face) exp = exp || exposed[i_face];
          if (exp) {
            bool bad = false;
            auto eval = [&](Mat<> direction) {
              Mat<> center = elem.tree->center() + elem.tree->nominal_size()/2*direction;
              Mat<> nearest = surf_geom->nearest_point(center, huge, 2*elem.tree->nominal_size()).point();
              double tol = .3;
              bad = bad || (nearest - center).normalized().dot(direction.normalized()) < -tol;
            };
            for (int i_face = 0; i_face < 2*nd; ++i_face) if (exposed[i_face]) {
              Mat<> direction = math::direction(nd, i_face).cast<double>();
              eval(direction);
              for (int j_face = 0; j_face < 2*(i_face/2); ++j_face) if (exposed[j_face]) {
                Mat<> dir = direction;
                dir(j_face/2) = math::sign(j_face%2);
                eval(dir);
                for (int k_face = 0; k_face < 2*(j_face/2); ++k_face) if (exposed[k_face]) {
                  dir(k_face/2) = math::sign(k_face%2);
                  eval(dir);
                }
              }
            }
            if (bad) {
              changed = true;
              #pragma omp atomic write
              elem.record = 2;
            }
          }
        }
      }
    }
  } while (changed);
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
        elem.tree.unpair();
      }
    }
  }
}

void Accessible_mesh::purge()
{
  if (tree) {
    // delete obsolete elements of `extrude_cons`
    erase_if(extrude_cons, [](Element_face_connection<Deformed_element>* con){return con->element(0).record == 2 || con->element(1).record == 2;});
    // delete connections to old elements (has to happen before deleting elements or else use after free)
    car.purge_connections();
    def.purge_connections();
    // delete old elements
    car.elems.purge();
    def.elems.purge();
  }
  // delete dangling vertex pointers
  erase_if(vert_ptrs, &Vertex::Non_transferable_ptr::is_null);
  for (auto& b_verts : boundary_verts) erase_if(b_verts, &Vertex::Non_transferable_ptr::is_null);
}

void Accessible_mesh::id_smooth_verts()
{
  smooth_verts.clear();
  auto verts = vertices();
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
    auto& vert = verts[i_vert];
    if (vert.is_mobile() && vert.needs_smooth()) smooth_verts.emplace_back(vert);
  }
}

bool Accessible_mesh::update(std::function<bool(Element&)> refine_criterion, std::function<bool(Element&)> unrefine_criterion)
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
    bool ref = refine_criterion(elem);
    bool unref = unrefine_criterion(elem);
    if (ref && !unref) elem.record = 1;
    else if (unref && !ref) elem.record = -1;
    elem.set_needs_smooth(false);
  }
  // pass refinement requests of extruded elements to their extrusion parents
  #pragma omp parallel for
  for (auto con : extrude_cons) {
    auto& inside = con->element(1);
    Lock::Acquire a(inside.lock);
    if (inside.record == 0) inside.record = con->element(0).record;
    else inside.record = std::max(inside.record, con->element(0).record);
    inside.unrefinement_locked = inside.unrefinement_locked || con->element(0).unrefinement_locked;
  }
  int n_orig [2];
  // refine elements
  for (bool is_deformed : {0, 1}) n_orig[is_deformed] = container(is_deformed).element_view().size(); // count how many elements there are before adding, so we know where the new ones start
  for (bool is_deformed : {0, 1}) refine_by_record(is_deformed, 0, n_orig[is_deformed]);
  bool changed;
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
            unref = !elem.unrefinement_locked;
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
  // incremental flood fill
  do {
    changed = false;
    // synchronize refinement level of surface elements with their non-surface neighbors
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      if (exists(elem.tree)) {
        for (int i_face = 0; i_face < 2*nd; ++i_face) {
          Tree* neighbor = elem.tree->find_neighbor(math::direction(nd, i_face));
          if (neighbor) {
            if (!exists(neighbor)) {
              if (neighbor->refinement_level() > elem.refinement_level()) {
                Tree* p = neighbor->parent();
                bool can_unref = true;
                for (Tree* child : p->children()) can_unref = can_unref && !exists(child);
                if (can_unref && !needs_refine(p)) {
                  changed = true;
                  for (Tree* child : p->children()) {
                    if (child->elem) {
                      child->elem->record = 2;
                    }
                  }
                  p->unrefine();
                }
              } else if (neighbor->refinement_level() < elem.refinement_level() - 1) {
                changed = true;
                refine_set_status(neighbor);
              } else if (neighbor->refinement_level() < elem.refinement_level()) {
                int min_rl = std::numeric_limits<int>::max();
                for (int j_face = 0; j_face < 2*nd; ++j_face) {
                  Tree* n = neighbor->find_neighbor(math::direction(nd, j_face));
                  if (exists(n)) min_rl = std::min(min_rl, n->refinement_level());
                }
                if (neighbor->refinement_level() < min_rl) {
                  changed = true;
                  refine_set_status(neighbor);
                }
              }
            }
          }
        }
      }
    }
    // add new elements
    for (bool is_deformed : {0, 1}) {
      auto& cont = container(is_deformed);
      auto& cont_elems = cont.element_view();
      int sz = cont_elems.size();
      for (int i_elem = 0; i_elem < sz; ++i_elem) {
        auto& elem = cont_elems[i_elem];
        if (exists(elem.tree)) {
          for (int i_face = 0; i_face < 2*nd; ++i_face) {
            for (Tree* neighbor : elem.tree->find_neighbors(math::direction(nd, i_face))) {
              if (!exists(neighbor)) if (!is_surface(neighbor)) {
                changed = true;
                add_elem(is_deformed, *neighbor).record = 0;
                neighbor->set_status(1);
              }
            }
          }
        }
      }
    }
  } while (changed);
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
  delete_bad_extrusions();
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
  int n_before = elems.size();
  purge();
  int n_after = elems.size();
  // connect new elements
  connect_new<         Element>(0);
  connect_new<Deformed_element>(0);
  extrude(true);
  connect_rest(surf_bc_sn);
  // if any immobile vertices have strayed from their nominal position (probably by `eat`ing)
  // snap them back where they belong
  auto& car_elems = cartesian().elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < car_elems.size(); ++i_elem) {
    auto& elem = car_elems[i_elem];
    for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
      Lock::Acquire a(elem.vertex(i_vert).lock);
      auto& pos = elem.vertex(i_vert).pos;
      double nom_sz = elem.nominal_size();
      auto nom_pos = elem.nominal_position();
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        pos[i_dim] = nom_sz*(nom_pos[i_dim] + (i_vert/math::pow(2, nd - 1 - i_dim))%2);
      }
      pos(Eigen::seqN(0, nd)) += elem.origin;
    }
  }
  id_boundary_verts();
  snap_vertices();
  id_smooth_verts();
  return n_before > n_after; // any change to the element structure (including adding elements!) will cause `purge` to reduce the size of `elems`
}

void Accessible_mesh::set_all_smooth()
{
  auto& elems = elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].set_needs_smooth(true);
  }
  id_smooth_verts();
  id_boundary_verts();
}

void Accessible_mesh::relax(double factor)
{
  id_boundary_verts();
  // calculate average neighbor position
  #pragma omp parallel for
  for (auto& vert : smooth_verts) {
    vert->temp_vector.setZero();
    auto neighbs = vert->get_neighbors();
    for (auto& neighb : neighbs) vert->temp_vector += neighb.pos;
    vert->temp_vector /= neighbs.size();
  }
  // update position
  #pragma omp parallel for
  for (auto& vert : smooth_verts) {
    if (vert->is_mobile()) vert->pos = factor*vert->temp_vector + (1 - factor)*vert->pos;
  }
  snap_vertices();
}

void Accessible_mesh::reset_verts()
{
  int nv = params.n_vertices();
  auto verts = vertices();
  #pragma omp parallel for
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
    verts[i_vert].temp_vector = verts[i_vert].pos;
  }
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements().size(); ++i_elem) {
    auto& elem = elements()[i_elem];
    if (elem.tree) {
      for (int i_vert = 0; i_vert < nv; ++i_vert) {
        auto& vert = elem.vertex(i_vert);
        Lock::Acquire a(vert.lock);
        vert.pos = elem.tree->nominal_position();
        for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
          vert.pos(i_dim) += elem.tree->nominal_size()*((i_vert/math::pow(2, params.n_dim - 1 - i_dim))%2);
        }
      }
    }
  }
  #pragma omp parallel for
  for (unsigned i_con = 0; i_con < extrude_cons.size(); ++i_con) {
    auto con = extrude_cons[i_con];
    auto& elem = con->element(0);
    auto dir = con->direction();
    int stride = math::pow(2, params.n_dim - 1 - dir.i_dim[0]);
    for (int i_vert = 0; i_vert < nv; ++i_vert) {
      auto& vert = elem.vertex(i_vert);
      Lock::Acquire a(vert.lock);
      int face_sign = (i_vert/stride)%2;
      if (face_sign == dir.face_sign[1]) {
        vert.pos = elem.vertex(i_vert + stride*(dir.face_sign[0] - face_sign)).pos;
      }
    }
  }
  verts_are_reset = true;
}

void Accessible_mesh::restore_verts()
{
  verts_are_reset = false;
  auto verts = vertices();
  #pragma omp parallel for
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
    verts[i_vert].pos = verts[i_vert].temp_vector;
  }
}

template <typename T>
void h5_write_row(H5::DataSet& dset, int cols, int i_row, T* data)
{
  hsize_t row_dims [2] {1, hsize_t(cols)};
  H5::DataSpace mspace (2, row_dims, nullptr);
  hsize_t offset [2] {hsize_t(i_row), 0};
  hsize_t stride [2] {1, 1};
  hsize_t block [2] {1, 1};
  auto dspace = dset.getSpace();
  dspace.selectHyperslab(H5S_SELECT_SET, row_dims, offset, stride, block);
  dset.write(data, dset.getDataType(), mspace, dspace);
}

template <typename T>
void h5_write_value(H5::DataSet& dset, int i_row, T data)
{
  h5_write_row(dset, 1, i_row, &data);
}

template <typename T>
void h5_add_attr(H5::H5Object& obj, std::string name, T value, H5::DataType dtype = H5::PredType::NATIVE_INT)
{
  hsize_t attr_dim = 1;
  H5::DataSpace dspace(1, &attr_dim);
  auto attr = obj.createAttribute(name, dtype, dspace);
  attr.write(dtype, &value);
}

template <typename T = int>
T h5_get_attr(H5::H5Object& obj, std::string name, H5::DataType dtype = H5::PredType::NATIVE_INT)
{
  auto attr = obj.openAttribute(name.c_str());
  T value;
  attr.read(dtype, &value);
  return value;
}

Storage_params read_params(std::string file_name)
{
  H5::H5File file(file_name + ".mesh.h5", H5F_ACC_RDONLY);
  Storage_params params {
    h5_get_attr(file, "n_stage"),
    h5_get_attr(file, "n_var"),
    h5_get_attr(file, "n_dim"),
    h5_get_attr(file, "row_size"),
    h5_get_attr(file, "n_forcing"),
  };
  return params;
}

double read_root_sz(std::string file_name)
{
  H5::H5File file(file_name + ".mesh.h5", H5F_ACC_RDONLY);
  return h5_get_attr<double>(file, "root_size", H5::PredType::NATIVE_DOUBLE);
}

template <typename T>
void h5_read_row(H5::DataSet& dset, int cols, int i_row, T* data)
{
  hsize_t row_dims [2] {1, hsize_t(cols)};
  H5::DataSpace mspace (2, row_dims, nullptr);
  hsize_t offset [2] {hsize_t(i_row), 0};
  hsize_t stride [2] {1, 1};
  hsize_t block [2] {1, 1};
  auto dspace = dset.getSpace();
  dspace.selectHyperslab(H5S_SELECT_SET, row_dims, offset, stride, block);
  dset.read(data, dset.getDataType(), mspace, dspace);
}

template <typename T>
T h5_read_value(H5::DataSet& dset, int i_row)
{
  T data;
  h5_read_row(dset, 1, i_row, &data);
  return data;
}

void Accessible_mesh::write(std::string name)
{
  H5::H5File file(name + ".mesh.h5", H5F_ACC_TRUNC);
  h5_add_attr(file, "version_major", config::version_major);
  h5_add_attr(file, "version_minor", config::version_minor);
  h5_add_attr(file, "version_patch", config::version_patch);
  h5_add_attr(file, "n_dim", params.n_dim);
  h5_add_attr(file, "row_size", params.row_size);
  h5_add_attr(file, "n_stage", params.n_stage);
  h5_add_attr(file, "n_var", params.n_var);
  h5_add_attr(file, "n_forcing", params.n_forcing);
  h5_add_attr(file, "root_size", root_sz, H5::PredType::NATIVE_DOUBLE);
  hsize_t dims[2];
  // write vertices
  file.createGroup("/vertices");
  auto verts = vertices();
  dims[0] = verts.size();
  dims[1] = 3;
  auto dset = file.createDataSet("vertices/position", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dims));
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
    auto& vert = verts[i_vert];
    h5_write_row(dset, 3, i_vert, vert.pos.data());
    vert.record.clear();
    vert.record.push_back(i_vert);
  }
  // write elements
  file.createGroup("/elements");
  dims[0] = elems.size();
  dims[1] = params.n_vertices();
  auto vert_dset = file.createDataSet("/elements/vertices", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
  dims[1] = params.n_dim;
  auto nom_pos_dset = file.createDataSet("/elements/nominal_position", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
  dims[1] = 1;
  auto is_def_dset = file.createDataSet("/elements/is_deformed", H5::PredType::NATIVE_HBOOL, H5::DataSpace(2, dims));
  auto ref_level_dset = file.createDataSet("/elements/refinement_level", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    int vert_inds [8] {};
    for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
      vert_inds[i_vert] = elem.vertex(i_vert).record[0];
    }
    h5_write_row(vert_dset, params.n_vertices(), i_elem, vert_inds);
    std::vector<int> nom_pos = elem.nominal_position();
    h5_write_row(nom_pos_dset, params.n_dim, i_elem, nom_pos.data());
    h5_write_value(is_def_dset, i_elem, elem.get_is_deformed());
    h5_write_value(ref_level_dset, i_elem, elem.refinement_level());
    elem.record = i_elem;
  }
  // write conformal connections
  dims[0] = car.cons.size() + def.cons.size();
  dims[1] = 6;
  file.createGroup("/connections");
  auto con_dset = file.createDataSet("/connections/conformal", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
  #define WRITE_CONS(start, cons) \
    for (int i_con = 0; i_con < int(cons.size()); ++i_con) { \
      auto& con = *cons[i_con]; \
      Connection_direction dir = con.get_direction(); \
      int data [6]; \
      for (int i_side = 0; i_side < 2; ++i_side) { \
        data[i_side] = con.element(i_side).record; /* record element indices */ \
        data[2 + i_side] = dir.i_dim[i_side]; \
        data[4 + i_side] = dir.face_sign[i_side]; \
      } \
      h5_write_row(con_dset, 6, start + i_con, data); \
    }
  WRITE_CONS(0, car.cons);
  WRITE_CONS(car.cons.size(), def.cons);
  #undef WRITE_CONS
  // write refined connections
  auto& car_cons = car.refined_connections();
  auto& def_cons = def.refined_connections();
  dims[0] = car_cons.size() + def_cons.size();
  dims[1] = 12;
  auto ref_con_dset = file.createDataSet("/connections/refined", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
  #define WRITE_REF_CONS(start, cons) \
    for (int i_con = 0; i_con < cons.size(); ++i_con) { \
      auto& con = cons[i_con]; \
      int data[12] {}; \
      data[0] = con.coarse_element().record; \
      for (int i_fine = 0; i_fine < con.n_fine_elements(); ++i_fine) { \
        data[1 + i_fine] = con.connection(i_fine).element(!con.order_reversed()).record; \
      } \
      data[5] = con.order_reversed(); \
      for (int i_dim : {0, 1}) data[6 + i_dim] = con.stretch()[i_dim]; \
      Con_dir<Deformed_element> dir = con.direction(); \
      for (int i_side = 0; i_side < 2; ++i_side) { \
        data[8 + i_side] = dir.i_dim[i_side]; \
        data[10 + i_side] = dir.face_sign[i_side]; \
      } \
      h5_write_row(ref_con_dset, 12, start + i_con, data); \
    }
  WRITE_REF_CONS(0, car_cons);
  WRITE_REF_CONS(car_cons.size(), def_cons);
  #undef WRITE_REF_CONS
  // write boundary connections
  auto& bound_cons = boundary_connections();
  dims[0] = bound_cons.size();
  dims[1] = 4;
  auto bound_con_dset = file.createDataSet("/connections/boundary", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
  for (int i_con = 0; i_con < bound_cons.size(); ++i_con) {
    auto& con = bound_cons[i_con];
    int data[4];
    data[0] = con.element().record;
    data[1] = con.bound_cond_serial_n();
    data[2] = con.i_dim();
    data[3] = con.inside_face_sign();
    h5_write_row(bound_con_dset, 4, i_con, data);
  }
  // write tree
  if (tree) {
    file.createGroup("/tree");
    dims[0] = 1;
    dims[1] = params.n_dim;
    auto orig_dset = file.createDataSet("/tree/origin", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dims));
    Mat<> origin = tree->origin();
    h5_write_row(orig_dset, params.n_dim, 0, origin.data());
    int n_vert = params.n_vertices();
    dims[0] = tree->count();
    dims[1] = 2 + n_vert;
    auto child_dset = file.createDataSet("/tree/children", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
    int row = 0;
    std::function<int(Tree*)> write_tree = [&](Tree* t) {
      std::vector<int> data(2 + n_vert, -1);
      if (t->elem) data[0] = t->elem->record;
      data[1] = t->get_status();
      int my_row = row++;
      auto children = t->children();
      for (unsigned i_child = 0; i_child < children.size(); ++i_child) data[2 + i_child] = write_tree(children[i_child]);
      h5_write_row(child_dset, 2 + n_vert, my_row, data.data());
      return my_row;
    };
    write_tree(tree.get());
  }
  // write face warping
  file.createGroup("/elements/face_warping");
  auto& def_elems = def.elements();
  dims[0] = def_elems.size();
  dims[1] = 1;
  auto ind_dset = file.createDataSet("/elements/face_warping/element_indices", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
  dims[1] = 2*params.n_dim*params.n_face_qpoint();
  auto adj_dset = file.createDataSet("/elements/face_warping/node_adjustments", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dims));
  for (int i_elem = 0; i_elem < def_elems.size(); ++i_elem) {
    auto& elem = def_elems[i_elem];
    h5_write_value(ind_dset, i_elem, elem.record);
    h5_write_row(adj_dset, dims[1], i_elem, elem.node_adjustments());
  }
}

void Accessible_mesh::read_file(std::string file_name)
{
  H5::H5File file(file_name + ".mesh.h5", H5F_ACC_RDONLY);
  hsize_t dims [2];
  // read elements
  auto vert_pos_dset = file.openDataSet("/vertices/position");
  auto vert_ind_dset = file.openDataSet("/elements/vertices");
  auto nom_pos_dset = file.openDataSet("/elements/nominal_position");
  auto is_def_dset = file.openDataSet("/elements/is_deformed");
  auto ref_level_dset = file.openDataSet("/elements/refinement_level");
  is_def_dset.getSpace().getSimpleExtentDims(dims);
  int n_elem = dims[0];
  int n_vert = params.n_vertices();
  std::vector<Element*> elem_ptrs(n_elem); // really need to switch to storing a flat array of fully polymorphic elements to avoid this nonsense
  std::vector<Deformed_element*> def_elem_ptrs(n_elem, nullptr);
  for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
    std::vector<int> nom_pos(params.n_dim);
    h5_read_row(nom_pos_dset, params.n_dim, i_elem, nom_pos.data());
    int ref_level = h5_read_value<int>(ref_level_dset, i_elem);
    int is_def = h5_read_value<bool>(is_def_dset, i_elem);
    int sn = add_element(ref_level, is_def, nom_pos, tree ? tree->origin() : Mat<>::Zero(params.n_dim));
    int vert_inds[8] {};
    h5_read_row(vert_ind_dset, n_vert, i_elem, vert_inds);
    auto& elem = element(ref_level, is_def, sn);
    elem_ptrs[i_elem] = &elem;
    if (is_def) def_elem_ptrs[i_elem] = &def.elems.at(ref_level, sn);
    for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
      h5_read_row(vert_pos_dset, params.n_dim, vert_inds[i_vert], elem.vertex(i_vert).pos.data());
    }
  }
  // read tree
  if (tree) {
    auto child_dset = file.openDataSet("/tree/children");
    std::function<void(Tree*, int)> read_tree = [&](Tree* t, int row) {
      std::vector<int> data(2 + n_vert);
      h5_read_row(child_dset, 2 + n_vert, row, data.data());
      if (data[0] >= 0) {
        Element& elem = *elem_ptrs[data[0]];
        t->elem.pair(elem.tree);
        t->def_elem = def_elem_ptrs[data[0]]; // if not deformed, this is just `nullptr`, as it should be
      }
      t->set_status(data[1]);
      if (data[2] >= 0) {
        t->refine();
        auto children = t->children();
        for (int i_child = 0; i_child < n_vert; ++i_child) read_tree(children[i_child], data[2 + i_child]);
      }
    };
    read_tree(tree.get(), 0);
  }
  // read conformal connections
  auto con_dset = file.openDataSet("/connections/conformal");
  con_dset.getSpace().getSimpleExtentDims(dims);
  int n_con = dims[0];
  for (int i_con = 0; i_con < n_con; ++i_con) {
    int data [6];
    h5_read_row(con_dset, 6, i_con, data);
    std::array<Element*, 2> el_ar;
    for (int i_side = 0; i_side < 2; ++i_side) el_ar[i_side] = elem_ptrs[data[i_side]];
    if (el_ar[0]->get_is_deformed() && el_ar[1]->get_is_deformed()) {
      std::array<Deformed_element*, 2> def_el_ar {def_elem_ptrs[data[0]], def_elem_ptrs[data[1]]};
      def.cons.emplace_back(new Element_face_connection<Deformed_element>(def_el_ar, {{data[2], data[3]}, {bool(data[4]), bool(data[5])}}));
      if (bool(def_el_ar[0]->tree) != bool(def_el_ar[1]->tree)) extrude_cons.push_back(def.cons.back().get());
    } else {
      car.cons.emplace_back(new Element_face_connection<Element>(el_ar, {data[2]}));
    }
  }
  // read refined connections
  auto ref_con_dset = file.openDataSet("/connections/refined");
  ref_con_dset.getSpace().getSimpleExtentDims(dims);
  n_con = dims[0];
  for (int i_con = 0; i_con < n_con; ++i_con) {
    int data [12];
    h5_read_row(ref_con_dset, 12, i_con, data);
    std::array<bool, 2> stretch {bool(data[6]), bool(data[7])};
    int n_fine = math::pow(2, params.n_dim - 1 - stretch[0] - stretch[1]);
    Element* coarse = elem_ptrs[data[0]];
    std::vector<Element*> fine(n_fine);
    bool is_def = coarse->get_is_deformed();
    for (int i_fine = 0; i_fine < n_fine; ++i_fine) {
      fine[i_fine] = elem_ptrs[data[1 + i_fine]];
      is_def = is_def && fine[i_fine]->get_is_deformed();
    }
    if (is_def) {
      Deformed_element* def_coarse = def_elem_ptrs[data[0]];
      std::vector<Deformed_element*> def_fine(n_fine);
      for (int i_fine = 0; i_fine < n_fine; ++i_fine) def_fine[i_fine] = def_elem_ptrs[data[1 + i_fine]];
      def.ref_face_cons[math::log(2, n_fine)].emplace_back(
        new Refined_connection<Deformed_element>{def_coarse, def_fine, {{data[8], data[9]}, {bool(data[10]), bool(data[11])}}, bool(data[5]), stretch}
      );
    } else {
      car.ref_face_cons[params.n_dim - 1].emplace_back(new Refined_connection<Element>{coarse, fine, {data[8]}, bool(data[5])});
    }
  }
  // read boundary connections
  auto bound_con_dset = file.openDataSet("/connections/boundary");
  bound_con_dset.getSpace().getSimpleExtentDims(dims);
  for (int i_con = 0; i_con < int(dims[0]); ++i_con) {
    int data [4];
    h5_read_row(bound_con_dset, 4, i_con, data);
    HEXED_ASSERT(data[1] < int(bound_conds.size()), "mesh file refers to nonexistant boundary condition");
    if (elem_ptrs[data[0]]->get_is_deformed()) {
      def.bound_cons.emplace_back(new Typed_bound_connection<Deformed_element>(*def_elem_ptrs[data[0]], data[2], data[3], data[1]));
    } else {
      car.bound_cons.emplace_back(new Typed_bound_connection<Element         >(    *elem_ptrs[data[0]], data[2], data[3], data[1]));
    }
  }
  // read face warping
  auto ind_dset = file.openDataSet("/elements/face_warping/element_indices");
  auto adj_dset = file.openDataSet("/elements/face_warping/node_adjustments");
  adj_dset.getSpace().getSimpleExtentDims(dims);
  for (unsigned row = 0; row < dims[0]; ++row) {
    Deformed_element* elem = def_elem_ptrs[h5_read_value<int>(ind_dset, row)];
    HEXED_ASSERT(elem, "file specifies face warping for a Cartesian element");
    h5_read_row(adj_dset, dims[1], row, elem->node_adjustments());
  }
  cleanup();
}

Accessible_mesh::Accessible_mesh(std::string file_name, std::vector<Flow_bc*> extremal_bcs, Surface_geom* geometry, Flow_bc* surface_bc)
: Accessible_mesh(read_params(file_name), read_root_sz(file_name))
{
  // take ownership of these to avoid memory leaks in case of exception
  std::unique_ptr<Flow_bc> fbc;
  if (surface_bc) fbc.reset(surface_bc);
  std::unique_ptr<Surface_geom> g;
  if (geometry) g.reset(geometry);
  // create the tree
  {
    H5::H5File file(file_name + ".mesh.h5", H5F_ACC_RDONLY);
    HEXED_ASSERT(file.exists("tree"), "attempt to read a non-tree mesh from a file as a tree mesh");
    Mat<> orig(params.n_dim);
    auto orig_dset = file.openDataSet("/tree/origin");
    h5_read_row(orig_dset, params.n_dim, 0, orig.data());
    create_tree(extremal_bcs, orig);
  }
  HEXED_ASSERT(bool(fbc) == bool(g), "must specify both surface geometry and surface boundary condition or neither");
  if (surface_bc) {
    surf_bc_sn = add_boundary_condition(fbc.release(), new Geom_mbc(g.release()));
    surf_geom = geometry;
  }
  read_file(file_name);
}

Accessible_mesh::Accessible_mesh(std::string file_name, std::vector<Flow_bc*> flow_bcs, std::vector<Mesh_bc*> mesh_bcs)
: Accessible_mesh(read_params(file_name), read_root_sz(file_name))
{
  HEXED_ASSERT(flow_bcs.size() == mesh_bcs.size(), "must supply same number of flow and mesh boundary conditions");
  for (unsigned i_bc = 0; i_bc < flow_bcs.size(); ++i_bc) {
    add_boundary_condition(flow_bcs[i_bc], mesh_bcs[i_bc]);
  }
  read_file(file_name);
}

void write_polymesh_file(std::string dir_name, std::string name, std::string cls, int n_entries, std::function<std::string(int)> entries, std::string note = "")
{
  std::ofstream file(dir_name + name);
  file
    << "// this file was generated for OpenFOAM by Hexed, an open-source mesher and CFD solver\n"
    << "// https://github.com/ARTLab-GT/hexed\n\n"
    << "FoamFile\n"
    << "{\n"
    << "    version 2.0;\n"
    << "    format ascii;\n";
  if (!note.empty()) file << "    note \"" << note << "\";\n";
  file
    << "    class " << cls << ";\n"
    << "    location \"constant/polyMesh\";\n"
    << "    object " << name << ";\n"
    << "}\n"
    << "\n" << n_entries << "\n(\n";
  for (int i_entry = 0; i_entry < n_entries; ++i_entry) file << "    " << entries(i_entry) << "\n";
  file << ")\n";
}

void Accessible_mesh::export_polymesh(std::string dir_name)
{
  dir_name = dir_name + "polyMesh/";
  if (std::filesystem::exists(dir_name)) std::filesystem::remove_all(dir_name);
  std::filesystem::create_directory(dir_name);
  auto verts = vertices();
  auto& bound_cons = boundary_connections();
  auto& elems = elements();
  auto& elem_cons = element_connections();
  int n_internal = elem_cons.size();
  int n_faces = n_internal + bound_cons.size();
  std::string face_note = format_str(200, "nPoints:%i nCells:%i nFaces:%i nInternalFaces:%i",
                                     verts.size(), elems.size(), n_faces, n_internal);
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) elems[i_elem].record = i_elem;
  write_polymesh_file(dir_name, "points", "vectorField", verts.size(), [&](int i_vert) {
    Vertex& vert = verts[i_vert];
    vert.record.clear();
    vert.record.push_back(i_vert);
    return format_str(100, "(%.20e %.20e %.20e)", vert.pos[0], vert.pos[1], vert.pos[2]);
  });
  std::vector<int> owners(n_faces);
  std::vector<int> neighbors(n_internal);
  std::vector<std::vector<int>> faces(n_faces);
  int i_face = 0;
  auto add_verts = [&](Element& elem, int i_dim, int face_sign, bool flip) {
    auto& verts = faces[i_face++];
    for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
      if ((i_vert/math::pow(2, params.n_dim - i_dim - 1))%2 == face_sign) {
        verts.push_back(elem.vertex(i_vert).record[0]);
      }
    }
    if ((face_sign != i_dim%2) != flip) std::swap(verts[0], verts[1]);
    else std::swap(verts[2], verts[3]);
  };
  for (int i_con = 0; i_con < elem_cons.size(); ++i_con) {
    auto& con = elem_cons[i_con];
    int owner_side = con.element(1).record < con.element(0).record;
    owners[i_face] = con.element(owner_side).record;
    neighbors[i_face] = con.element(!owner_side).record;
    int finer_side = con.element(1).refinement_level() > con.element(0).refinement_level();
    auto dir = con.get_direction();
    add_verts(con.element(finer_side), dir.i_dim[finer_side], dir.face_sign[finer_side], finer_side != owner_side);
  }
  std::vector<int> n_bound_cons(bound_conds.size(), 0);
  std::vector<int> bc_starts(bound_conds.size(), n_internal);
  std::vector<std::string> bc_names(bound_conds.size());
  std::vector<std::string> bc_types(bound_conds.size());
  for (int i_bc = 0; i_bc < int(bound_conds.size()); ++i_bc) {
    for (int i_con = 0; i_con < bound_cons.size(); ++i_con) {
      auto& con = bound_cons[i_con];
      if (con.bound_cond_serial_n() == i_bc) {
        ++n_bound_cons[i_bc];
        owners[i_face] = con.element().record;
        add_verts(con.element(), con.i_dim(), con.inside_face_sign(), false);
      }
    }
    if (i_bc) bc_starts[i_bc] = bc_starts[i_bc - 1] + n_bound_cons[i_bc - 1];
    if (tree) {
      if (i_bc == surf_bc_sn) {
        bc_names[i_bc] = "surface_bc";
        bc_types[i_bc] = "wall";
      } else {
        bc_names[i_bc] = format_str(100, "extremal_bc%i%i", i_bc/2, i_bc%2);
        bc_types[i_bc] = "patch";
      }
    } else {
      bc_names[i_bc] = format_str(100, "bc%i", i_bc);
      bc_types[i_bc] = "patch";
    }
  }
  write_polymesh_file(dir_name, "faces", "faceList", n_faces, [&](int i_entry) {
    auto& verts = faces[i_entry];
    return format_str(100, "4(%i %i %i %i)", verts[0], verts[1], verts[2], verts[3]);
  });
  write_polymesh_file(dir_name, "boundary", "polyBoundaryMesh", bound_conds.size(), [&](int i_bc) {
    return format_str(200, "%s {type %s; nFaces %i; startFace %i;}", bc_names[i_bc].c_str(), bc_types[i_bc].c_str(), n_bound_cons[i_bc], bc_starts[i_bc]);
  });
  write_polymesh_file(dir_name, "owner",     "labelList", n_faces,    [&](int i_entry){return format_str(100, "%i", owners   [i_entry]);}, face_note);
  write_polymesh_file(dir_name, "neighbour", "labelList", n_internal, [&](int i_entry){return format_str(100, "%i", neighbors[i_entry]);}, face_note);
}

}
