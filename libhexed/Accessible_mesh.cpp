#include <filesystem>
#include <fstream>
#include <queue>
#include <H5Cpp.h>
#include <hexed/Accessible_mesh.hpp>
#include <hexed/math.hpp>
#include <hexed/Row_index.hpp>
#include <hexed/erase_if.hpp>
#include <hexed/utils.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/History_monitor.hpp>
#include <hexed/Printer.hpp>
#include <hexed/Visualizer.hpp>

namespace hexed {

// makes an assertion and if it fails, visualizes the mesh before throwing
#define VIS_ASSERT(expression, message, ...) \
    HEXED_ASSERT(expression, \
                 _vis_return(message + std::string(" Writing diagnostic visualization to `meshing_diagnostic.*`")) \
                 __VA_OPT__(,) __VA_ARGS__) \

std::string Accessible_mesh::_vis_return(std::string str) {
  auto blocks = _blocks.boundary_sides();
  #pragma omp parallel for
  for (auto& b : blocks) b.reset();
  visualize("default", "meshing_diagnostic");
  return str;
}

Element_container& Accessible_mesh::container(bool is_deformed) {
  Element_container* containers [] {&car.elems, &def.elems};
  return *containers[is_deformed];
}

template<> Mesh_by_type<         Element>& Accessible_mesh::mbt() {return car;}
template<> Mesh_by_type<Deformed_element>& Accessible_mesh::mbt() {return def;}

void Accessible_mesh::_record_connections() {
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
}

namespace dijkstra {
  struct Node {
    next::Vertex* vert;
    double cost;
    int updates;
  };
  bool compare(Node x, Node y) {
    return x.cost > y.cost;
  }
}

void Accessible_mesh::_offset_vertices(double offset) {
  int nd = params.n_dim;
  auto verts = _blocks.verts();
  #pragma omp parallel for
  for (auto& vert : verts) {
    vert.offset.setZero();
    vert.dijkstra_dist = 0;
  }
  int nv = params.n_vertices()/2;
  for (auto& con : def.cons) {
    auto dir = con->get_direction();
    if (dir.i_dim[0] != dir.i_dim[1]) continue;
    bool is_new [2] {false, false};
    for (int i_side = 0; i_side < 2; ++i_side) {
      is_new[i_side] = con->element(i_side).active_shape().is_new;
    }
    if (is_new[0] != is_new[1]) {
      int new_elem = is_new[1];
      auto i_verts = vertex_inds(nd, dir)[0];
      Mat<3, dyn> vert_pos(3, nv);
      std::vector<next::Vertex*> con_verts(nv);
      for (int i_vert = 0; i_vert < nv; ++i_vert) {
        con_verts[i_vert] = &con->element(0).active_shape().vertex(i_verts[i_vert]);
        vert_pos(all, i_vert) = con_verts[i_vert]->unwarped_point();
      }
      for (int i_vert = 0; i_vert < nv; ++i_vert) {
        Mat<3, 2> edges;
        edges(all, 1).setUnit(2);
        for (int i_dim = 0; i_dim < nd - 1; ++i_dim) {
          int stride = math::pow(2, nd - 2 - i_dim);
          int start = i_vert - i_vert/stride%2*stride;
          edges(all, i_dim) = vert_pos(all, start + stride) - vert_pos(all, start);
        }
        Mat<3> nrml = edges(all, 0).cross(edges(all, 1)).normalized()
                      *math::sign(dir.face_sign[0])*math::sign(new_elem)*math::sign(dir.i_dim[0] == 1);
        double dot = nrml.dot(con_verts[i_vert]->offset);
        Mat<3> diff = nrml;
        if (dot < 0) {
          double norm_sq = con_verts[i_vert]->offset.squaredNorm();
          if (norm_sq > .1) { // `offset` should be 0 or >= 1
            diff -= con_verts[i_vert]->offset*dot/norm_sq;
            diff /= diff.dot(nrml);
          }
        }
        con_verts[i_vert]->offset += std::max(0., 1 - dot)*diff;
      }
    }
  }
  #pragma omp parallel for
  for (auto& vert : verts) {
    vert.set_pos(vert.unwarped_point() + offset*vert.nominal_size()*vert.offset);
  }
}

Mat<3> Accessible_mesh::_get_snapping_target(next::Vertex& vert, Mat<3> pos) {
  HEXED_ASSERT(Int(vert.record.size()) == 2*params.n_dim + 1, "Vertex record has not been set correctly.");
  auto seq = Eigen::seqN(0, params.n_dim);
  double ns = vert.nominal_size();
  if (vert.record[2*params.n_dim]) {
    if (vert.snapped_point >= 0) {
      return surf_geom->points()[vert.snapped_point];
    } else if (vert.snapped_edge >= 0) {
      auto& geom_edge = surf_geom->edges()[vert.snapped_edge];
      Array<double> nodes{geom_edge.nodes()};
      Int n_points = nodes.shape()[0];
      if (vert.snapped_endpoint == -1) {
        Int nearest = geom_edge.nearest_point(pos(seq), huge).index;
        if (nearest >= 0) pos = nodes(nearest).vector();
      } else {
        pos = nodes(vert.snapped_endpoint*(n_points - 1)).vector();
      }
    } else {
      #if 1
      pos(seq) = surf_geom->nearest_point(pos(seq), huge, ns/2).point();
      #else
      Mat<3> p0 = pos;
      bool found = false;
      for (auto n : vert.neighbors()) if (n) {
        if (!n->is_surface()) {
          p0 = n->unwarped_point();
          found = true;
        }
      }
      HEXED_ASSERT(found, "no non-surface neighbor found")
      auto sects = surf_geom->intersections(p0(seq), pos(seq));
      double best_sect = huge;
      found = false;
      for (double sect : sects) {
        if (sect > 0 && sect < best_sect) {
          best_sect = sect;
          found = true;
        }
      }
      if (found) pos = best_sect*pos + (1 - best_sect)*p0;
      #endif
    }
  }
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    for (int sign : {0, 1}) {
      if (vert.record[2*i_dim + sign]) {
        pos(i_dim) = tree->origin()(i_dim) + sign*tree->nominal_size();
      }
    }
  }
  return pos;
}

void Accessible_mesh::_fit_surface() {
  if (!surf_geom) return;
  Stopwatch_tree::Starter sw_fit(_stopwatch["update"]["fit surface"]);
  _blocks.edges_2d();
  _blocks.faces_3d();
  auto all_verts = _blocks.verts();
  auto& elems = def.elements();
  #pragma omp parallel for
  for (Int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    if (!elems[i_elem].tree) elems[i_elem].active_shape().is_new = true;
  }
  #pragma omp parallel for
  for (auto& vert : all_verts) {
    vert.set_pos(vert.nominal_position());
  }
  _offset_vertices(.2);
  #if 1
  {
    Task_message message(printers::info, "  Pre-edge-matching mesh optimization", "\n", "  ");
    _optimize(1, 10, true);
  }
  #endif
  #pragma omp parallel for
  for (auto& vert : all_verts) {
    vert.record.clear();
  }
  auto verts = _blocks.boundary_verts();
  #pragma omp parallel for
  for (next::Vertex& vert : all_verts) {
    Mat<3> point = vert.unwarped_point();
    vert.dijkstra_point = point;
    vert.snapped_point = -1;
    vert.snapped_edge = -1;
    vert.snapped_endpoint = -1;
  }
  auto faces = _blocks.faces_3d();
  #pragma omp parallel for
  for (next::Face& face : faces) {
    for (int i_edge = 0; i_edge < 4; ++i_edge) {
      face.edge(i_edge).snapped_edge = -1;
    }
  }
  auto find_nearest_vert = [&](Mat<3> point, int i_geom_edge = -1)->next::Vertex* {
    next::Vertex* nearest_vert = nullptr;
    double dist_sq = huge;
    for (auto& vert : verts) {
      if (!vert.glued()) {
        double ns = vert.nominal_size();
        double d = (vert.dijkstra_point - point).squaredNorm();
        for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
          for (int sign : {0, 1}) {
            double extreme = tree->origin()(i_dim) + sign*tree->nominal_size();
            if ((std::abs(point(i_dim) - extreme) > 3e-2*ns) !=
                (std::abs(vert.dijkstra_point(i_dim) - extreme) > 3e-2*ns)) d = huge;
          }
        }
        // don't bother to account for snapped neighbors unless d is initially < dist_sq
        if (d < std::min(ns, dist_sq)) {
          bool snapped_neighbor = false;
          for (next::Vertex* v : vert.neighbors()) {
            if (v) snapped_neighbor = snapped_neighbor || (v->snapped_edge != -1 && v->snapped_edge != i_geom_edge)
                                                       || v->snapped_point != -1;
          }
          if (snapped_neighbor) d *= 10;
          if (d < dist_sq) { // now we know d accounting for snapped neighbors, so this is the real comparison
            dist_sq = d;
            nearest_vert = &vert;
          }
        }
      }
    }
    return nearest_vert;
  };

  auto edges = surf_geom->edges();
  if (params.n_dim == 3) {
    for (Int i_geom_edge = 0; i_geom_edge < (Int)edges.size(); ++i_geom_edge) {
      auto& geom_edge = edges[i_geom_edge];
      std::array<next::Vertex*, 2> start_end {nullptr, nullptr};
      for (int i_endpoint = 0; i_endpoint < 2; ++i_endpoint) {
        Mat<3> endpoint = geom_edge.nodes()(i_endpoint*(geom_edge.nodes().shape()[0] - 1)).vector();
        start_end[i_endpoint] = find_nearest_vert(endpoint, i_geom_edge);
        if (start_end[i_endpoint]) {
          start_end[i_endpoint]->dijkstra_point = endpoint;
          start_end[i_endpoint]->snapped_edge = i_geom_edge;
          start_end[i_endpoint]->snapped_endpoint = i_endpoint;
        }
      }
      if (!start_end[0] || !start_end[1] || start_end[0] == start_end[1]) continue;
      #pragma omp parallel for
      for (next::Vertex& vert : verts) {
        vert.dijkstra_dist = huge;
        vert.dijkstra_updates = 0;
        vert.dijkstra_prev_vert = nullptr;
        vert.dijkstra_prev_edge = nullptr;
        double d = huge;
        auto nearest = geom_edge.nearest_point(vert.dijkstra_point, d);
        if (nearest.index >= 0 && nearest.distance <= d) {
          vert.dijkstra_curve_dist_sq = nearest.distance*nearest.distance;
          if (vert.snapped_edge >= 0) {
            vert.dijkstra_curve_dist_sq *= 100;
          }
          vert.dijkstra_arc_len = geom_edge.arc_length()[nearest.index];
        } else {
          vert.dijkstra_curve_dist_sq = std::sqrt(huge);
          vert.dijkstra_arc_len = std::sqrt(huge);
        }
      }
      std::priority_queue<
        dijkstra::Node,
        std::vector<dijkstra::Node>,
        std::function<bool(dijkstra::Node, dijkstra::Node)>
      > unvisited(&dijkstra::compare);
      // Dijkstra's algorithm will start at the second endpoint and go to the first,
      // so that we can traverse the path in reverse via `Vertex::dijkstra_prev`,
      // we will end up with a path from the first endpoint to the second
      unvisited.emplace(start_end[1], 0., 1);
      start_end[1]->dijkstra_dist = 0;
      start_end[1]->dijkstra_updates = 1;
      dijkstra::Node curr {nullptr, 0., 0};
      while (curr.vert != start_end[0] && !unvisited.empty()) {
        curr = unvisited.top();
        unvisited.pop();
        if (curr.updates < curr.vert->dijkstra_updates) continue;
        for (next::Edge& edge : curr.vert->edges()) if (!edge.glued()) {
          next::Vertex* vert = &edge.vertex(&edge.vertex(0) == curr.vert);
          double interval = std::max((curr.vert->dijkstra_point - vert->dijkstra_point).norm(),
                                     std::abs(vert->dijkstra_arc_len - curr.vert->dijkstra_arc_len));
          double d = curr.cost + .5*(curr.vert->dijkstra_curve_dist_sq + vert->dijkstra_curve_dist_sq)*interval;
          if (d < vert->dijkstra_dist) {
            vert->dijkstra_dist = d;
            vert->dijkstra_prev_vert = curr.vert;
            vert->dijkstra_prev_edge = &edge;
            ++vert->dijkstra_updates;
            unvisited.emplace(vert, d, vert->dijkstra_updates);
          }
        }
      }
      matched_vertices[i_geom_edge].clear();
      matched_edges[i_geom_edge].clear();
      if (curr.vert == start_end[0]) {
        next::Vertex* vert = curr.vert;
        do {
          matched_vertices[i_geom_edge].emplace_back(vert);
          if (vert->dijkstra_prev_edge) {
            vert->dijkstra_prev_edge->snapped_edge = i_geom_edge;
          }
          if (vert->snapped_edge < 0) {
            Int ind = geom_edge.nearest_point(vert->dijkstra_point, huge).index;
            if (ind >= 0) vert->dijkstra_point = geom_edge.nodes()(ind).vector();
            vert->snapped_edge = i_geom_edge;
          }
          vert = vert->dijkstra_prev_vert;
        } while (vert);
      }
    }
  } else if (params.n_dim == 2) {
    auto points = surf_geom->points();
    for (int i_point = 0; i_point < points.size(); ++i_point) {
      Mat<3> point {points[i_point][0], points[i_point][1], 0.};
      next::Vertex* vert = find_nearest_vert(point);
      if (vert) vert->snapped_point = i_point;
    }
  }
  for (auto& vert : verts) {
    if (vert.snapped_endpoint >= 0) {
      int n_snapped_edges = 0;
      for (auto& edge : vert.edges()) {
        n_snapped_edges += edge.snapped_edge != -1;
      }
      if (n_snapped_edges == 1) {
        for (auto& elem : vert.elements()) {
          auto face = elem.boundary_face_3d();
          if (face) {
            for (int i_edge = 0; i_edge < 4; ++i_edge) {
              auto* edge = &face->edge(i_edge);
              if (edge->glued()) edge = edge->glued_to();
              if (edge->snapped_edge == -1) edge->snapped_edge = -2;
              for (int i_vert = 0; i_vert < 2; ++i_vert) {
                if (edge->vertex(i_vert).snapped_edge == -1) edge->vertex(i_vert).snapped_edge = -2;
              }
            }
            break;
          }
        }
      }
    }
  }

  #pragma omp parallel for
  for (Int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    if (!elems[i_elem].tree) elems[i_elem].active_shape().is_new = true;
  }
  #pragma omp parallel for
  for (auto& vert : all_verts) {
    vert.set_pos(vert.nominal_position());
  }
  _offset_vertices(.2);
  #pragma omp parallel for
  for (Int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].active_shape().is_new = false;
  }
  Int elems_sz = elems.size();
  #pragma omp parallel for
  for (Int i_elem = 0; i_elem < elems_sz; ++i_elem) elems[i_elem].record = 0;

  if (params.n_dim == 3) {
    for (Int i_element = 0; i_element < elems_sz; ++i_element) {
      auto& elem = elems[i_element];
      for (int i_face = 0; i_face < 6; ++i_face) elem.face_record[i_face] = -1;
      const next::Element_shape* shape = elem.fake_shape();
      if (shape) {
        const next::Face* face = shape->boundary_face_3d();
        if (face) {
          int bf = shape->boundary_face();
          int i_dim = bf/2;
          bool i_sign = bf%2;
          std::vector<Int> matched_to(4);
          bool matched = false;
          for (int i_vert = 0; i_vert < 8; ++i_vert) if (i_vert/math::pow(2, 2 - i_dim)%2 == i_sign) {
            matched = matched || shape->vertex(i_vert).snapped_edge != -1;
          }
          for (int i_edge = 0; i_edge < 4; ++i_edge) {
            if (face->edge(i_edge).glued()) {
              for (int i_vert = 0; i_vert < 2; ++i_vert) {
                matched = matched || face->edge(i_edge).glued_to()->vertex(i_vert).snapped_edge != -1;
              }
            }
          }
          if (!matched) continue;
          for (int i_edge = 0; i_edge < 4; ++i_edge) {
            auto& edge = face->edge(i_edge);
            if (edge.glued()) {
              matched_to[i_edge] = edge.glued_to()->snapped_edge;
            } else {
              matched_to[i_edge] = edge.snapped_edge;
            }
          }
          auto set_vertices = [&](Element& e) {
            auto& s = e.shape();
            for (int i_vert = 0; i_vert < 8; ++i_vert) {
              s.vertex(i_vert).set_pos(shape->vertex(i_vert).unwarped_point());
            }
            for (int i_face = 0; i_face < 6; ++i_face) e.face_record[i_face] = -1;
            s.extruded_direction = shape->extruded_direction;
          };
          Int inside_sn = add_element(elem.refinement_level(), true, elem.nominal_position(), tree->origin(), 0);
          Deformed_element& inside = def.elems.at(elem.refinement_level(), inside_sn);
          set_vertices(inside);
          inside.shape().is_new = false;
          elem.face_record[2*i_dim + !i_sign] = inside_sn;
          Int surface_sn = add_element(elem.refinement_level(), true, elem.nominal_position(), tree->origin(), 0, bf);
          Deformed_element& surface = def.elems.at(elem.refinement_level(), surface_sn);
          set_vertices(surface);
          surface.shape().is_new = false;
          _connect({&inside, &surface}, Con_dir<Deformed_element>({i_dim, i_dim}, {i_sign, !i_sign}));
          elem.record = 2;
          std::vector<Deformed_element*> matched_elems(6, nullptr);
          for (int i_vert = 0; i_vert < 8; ++i_vert) {
            HEXED_ASSERT(std::isfinite(inside.shape().vertex(i_vert).unwarped_point().squaredNorm()),
                         "Vertex pos is not finite.")
            HEXED_ASSERT(std::isfinite(surface.shape().vertex(i_vert).unwarped_point().squaredNorm()),
                         "Vertex pos is not finite.")
          }
          for (int j_dim = 0; j_dim < 3; ++j_dim) if (j_dim != i_dim) {
            for (bool j_sign : {0, 1}) {
              int k_dim = 3 - j_dim - i_dim;
              int i_edge_matched = 2*(j_dim > k_dim) + j_sign;
              Int m = matched_to[i_edge_matched];
              if (m != -1) {
                Int sn = add_element(elem.refinement_level(), true, elem.nominal_position(), tree->origin(), 0, bf);
                Deformed_element& match_elem = def.elems.at(elem.refinement_level(), sn);
                set_vertices(match_elem);
                match_elem.shape().is_new = true;
                _connect({&surface, &match_elem}, Con_dir<Deformed_element>({j_dim, j_dim}, {j_sign, !j_sign}));
                _connect({&inside,  &match_elem}, Con_dir<Deformed_element>({j_dim, i_dim}, {j_sign, !i_sign}));
                matched_elems[2*j_dim + j_sign] = &match_elem;
                elem.face_record[2*j_dim + j_sign] = sn;
                for (bool k_sign : {0, 1}) {
                  int i_vert =   i_sign*math::pow(2, 2 - i_dim)
                               + j_sign*math::pow(2, 2 - j_dim)
                               + k_sign*math::pow(2, 2 - k_dim);
                  int i_snapped = shape->vertex(i_vert).snapped_edge;
                  HEXED_ASSERT(i_snapped != -1 || face->edge(i_edge_matched).glued(),
                               "vertex and edge do not agree on whether they are snapped")
                  auto& vert = match_elem.shape().vertex(i_vert);
                  vert.snapped_edge = i_snapped;
                  vert.snapped_endpoint = shape->vertex(i_vert).snapped_endpoint;
                  if (i_snapped >= 0) matched_vertices[i_snapped].emplace_back(&vert);
                }
                auto& matched_edge = match_elem.shape().boundary_face_3d()->edge(i_edge_matched);
                matched_edge.snapped_edge = m;
                if (m >= 0) matched_edges[m].emplace_back(&matched_edge);
                for (int i_vert = 0; i_vert < 8; ++i_vert) {
                  HEXED_ASSERT(std::isfinite(match_elem.shape().vertex(i_vert).unwarped_point().squaredNorm()),
                               "Vertex pos is not finite.")
                }
              } else {
                elem.face_record[2*j_dim + j_sign] = inside_sn;
                inside.face_record[2*j_dim + j_sign] = surface_sn;
              }
            }
          }
          Mat<3, 8> orig_pos;
          for (int i_vert = 0; i_vert < 8; ++i_vert) {
            orig_pos(all, i_vert) = shape->vertex(i_vert).unwarped_point();
          }
          for (int i_vert = 0; i_vert < 8; ++i_vert) {
            Mat<3> pos = orig_pos(all, i_vert);
            for (int j_dim = 0; j_dim < 3; ++j_dim) {
              double offset = .0 + .5*(j_dim == i_dim);
              int stride = math::pow(2, 2 - j_dim);
              bool j_sign = i_vert/stride%2;
              if ((j_dim == i_dim && j_sign != i_sign) || matched_elems[2*j_dim + j_sign]) {
                pos += offset*(orig_pos(all, i_vert - math::sign(j_sign)*stride) - orig_pos(all, i_vert));
              }
            }
            surface.shape().vertex(i_vert).set_pos(pos);
          }
          for (int j_dim = 0; j_dim < 3; ++j_dim) if (j_dim != i_dim) {
            int k_dim = 3 - j_dim - i_dim;
            for (bool j_sign : {0, 1}) if (matched_elems[2*j_dim + j_sign]) {
              for (bool k_sign : {0, 1}) {
                int i_vert =   i_sign*math::pow(2, 2 - i_dim)
                             + j_sign*math::pow(2, 2 - j_dim)
                             + k_sign*math::pow(2, 2 - k_dim);
                std::vector<Int> record {elem.refinement_level(), elem.face_record[2*j_dim + j_sign], k_dim, k_sign, i_dim, i_sign};
                auto& vert = surface.shape().vertex(i_vert);
                vert.record.insert(vert.record.end(), record.begin(), record.end());
              }
            }
          }
          std::vector<Deformed_element*> check_elems {&elem, &inside, &surface};
          for (auto p : matched_elems) if (p) check_elems.push_back(p);
          for (auto e : check_elems) {
            for (int i_vert = 0; i_vert < 8; ++i_vert) {
              HEXED_ASSERT(std::isfinite(e->shape().vertex(i_vert).unwarped_point().squaredNorm()),
                           "Vertex pos is not finite.")
            }
          }
        }
      }
    }
    for (int i_element = 0; i_element < elems.size(); ++i_element) {
      if (elems[i_element].record == 2) continue;
      for (int i_vert = 0; i_vert < 8; ++i_vert) {
        HEXED_ASSERT(std::isfinite(elems[i_element].shape().vertex(i_vert).unwarped_point().squaredNorm()),
                     "Vertex pos is not finite.")
        if (elems[i_element].fake_shape()) {
          HEXED_ASSERT(std::isfinite(elems[i_element].fake_shape()->vertex(i_vert).unwarped_point().squaredNorm()),
                       "Vertex pos is not finite.")
        }
      }
    }
    for (Int i_element = 0; i_element < elems.size(); ++i_element) {
      if (elems[i_element].record == 2) {
        elems[i_element].destroy_shape();
      } else {
        for (int i_vert = 0; i_vert < 8; ++i_vert) {
          Mat<3> p = elems[i_element].active_shape().vertex(i_vert).unwarped_point();
          elems[i_element].shape().vertex(i_vert).set_pos(p);
        }
      }
    }
    extrude_cons.clear();
    _blocks.boundary_sides();
    Int cons_sz = def.cons.size();
    Int ref_cons_sz = def.ref_face_cons[1].size();
    Int bound_cons_sz = def.bound_cons.size();
    for (Int i_con = 0; i_con < cons_sz; ++i_con) {
      auto& con = def.cons[i_con];
      if (!con) continue;
      auto dir = con->direction();
      bool replace = false;
      std::array<Deformed_element*, 2> elem_arr;
      std::array<Deformed_element*, 2> surfaces {nullptr, nullptr};
      for (int i_side = 0; i_side < 2; ++i_side) {
        Deformed_element& elem = con->element(i_side);
        if (elem.face_record[dir.i_face(i_side)] >= 0) {
          replace = true;
          elem_arr[i_side] = &def.elems.at(elem.refinement_level(), elem.face_record[dir.i_face(i_side)]);
          if (elem_arr[i_side]->face_record[dir.i_face(i_side)] >= 0) {
            surfaces[i_side] = &def.elems.at(elem.refinement_level(), elem_arr[i_side]->face_record[dir.i_face(i_side)]);
          }
        } else {
          elem_arr[i_side] = &elem;
        }
      }
      if (replace) {
        con.reset();
        if (bool(elem_arr[0]->fake_shape()) == bool(elem_arr[1]->fake_shape())) {
          _connect(elem_arr, dir);
          if (surfaces[0]) {
            HEXED_ASSERT(surfaces[1], "Both faces must identify a surface element or neither.");
            _connect(surfaces, dir);
          }
        } else {
          bool coarse_sign = surfaces[0];
          std::vector<Deformed_element*> fine {elem_arr[!coarse_sign], surfaces[!coarse_sign]};
          HEXED_ASSERT(fine[1], "Fine element must identify a surface element");
          std::array<bool, 2> stretch {false, false};
          HEXED_ASSERT(elem_arr[coarse_sign]->fake_shape(), "Coarse element must have a fake shape.");
          int i_bf = fine[1]->shape().boundary_face();
          int i_dim = i_bf/2;
          stretch[i_dim < 3 - i_dim - dir.i_dim[!coarse_sign]] = true;
          if (!(i_bf%2)) std::swap(fine[0], fine[1]);
          Con_dir<Deformed_element> new_dir {{dir.i_dim[coarse_sign], dir.i_dim[!coarse_sign]},
                                             {dir.face_sign[coarse_sign], dir.face_sign[!coarse_sign]}};
          _connect(elem_arr[coarse_sign], fine, new_dir, stretch);
        }
      }
    }
    for (Int i_con = 0; i_con < ref_cons_sz; ++i_con) {
      auto& con = def.ref_face_cons[1][i_con];
      if (!con) continue;
      auto dir = con->direction();
      bool reverse = con->order_reversed();
      Deformed_element* coarse;
      coarse = &con->coarse_element();
      bool replace = false;
      Int rec = coarse->face_record[dir.i_face(reverse)];
      if (rec >= 0) {
        coarse = &def.elems.at(coarse->refinement_level(), rec);
        replace = true;
      }
      std::vector<next::Element_shape*> coarse_shapes(4, &coarse->active_shape());
      rec = coarse->face_record[dir.i_face(reverse)];
      [[maybe_unused]] bool coarse_surface = false;
      if (rec >= 0) {
        Deformed_element* surface = &def.elems.at(coarse->refinement_level(), rec);
        int bf = surface->active_shape().boundary_face();
        for (int i_elem = 0; i_elem < 4; ++i_elem) {
          if (math::row_coordinate(2, 2, bf/2 > 3 - bf/2 - dir.i_dim[reverse], i_elem) == bf%2) {
            coarse_shapes[i_elem] = &surface->active_shape();
            coarse_surface = true;
          }
        }
      }
      auto stretch = con->stretch();
      std::vector<next::Element_shape*> fine_shapes;
      [[maybe_unused]] bool fine_surface = false;
      for (int i = 0, i_elem = 0; i < 1 + stretch[0]; ++i) {
        for (int i_fine = 0; i_fine < con->n_fine_elements(); ++i_fine) {
          Deformed_element* fine = &con->connection(i_fine).element(!reverse);
          rec = fine->face_record[dir.i_face(!reverse)];
          if (rec >= 0) {
            fine = &def.elems.at(fine->refinement_level(), rec);
            replace = true;
          }
          for (int j = 0; j < 1 + stretch[1]; ++j) {
            Deformed_element* f = fine;
            rec = f->face_record[dir.i_face(!reverse)];
            if (rec >= 0) {
              Deformed_element* surface = &def.elems.at(fine->refinement_level(), rec);
              int bf = surface->active_shape().boundary_face();
              if (bf >= 0) if (math::row_coordinate(2, 2, bf/2 > 3 - bf/2 - dir.i_dim[!reverse], i_elem) == bf%2) {
                f = surface;
                fine_surface = true;
              }
            }
            fine_shapes.push_back(&f->active_shape());
            ++i_elem;
          }
        }
      }
      if (replace) {
        //HEXED_ASSERT(fine_surface || !coarse_surface, "2-on-2 connection", assert::Not_implemented_error);
        Con_dir<Deformed_element> new_dir {{dir.i_dim[reverse], dir.i_dim[!reverse]},
                                           {dir.face_sign[reverse], dir.face_sign[!reverse]}};
        con.reset();
        next::Element_shape::connect({coarse_shapes, fine_shapes}, new_dir);
      }
    }
    for (Int i_con = 0; i_con < bound_cons_sz; ++i_con) {
      auto& con = def.bound_cons[i_con];
      if (!con) continue;
      auto dir = con->direction();
      Int record = con->element().face_record[dir.i_face(0)];
      if (record >= 0) {
        int bc_sn = con->bound_cond_serial_n();
        int ref_level = con->element().refinement_level();
        con.reset();
        connect_boundary(ref_level, true, record, dir.i_dim[0], dir.face_sign[0], bc_sn);
      }
    }
    for (auto& vert : all_verts) {
      if (vert.record.size() == 12) {
        std::array<Deformed_element*, 2> elem_arr;
        std::array<int, 2> dim_arr;
        std::array<bool, 2> sign_arr;
        for (int i_side = 0; i_side < 2; ++i_side) {
          elem_arr[i_side] = &def.elems.at(vert.record[6*i_side], vert.record[6*i_side + 1]);
          dim_arr[i_side] = vert.record[6*i_side + 2];
          sign_arr[i_side] = vert.record[6*i_side + 3];
        }
        int rotate = 0;
        for (int r : {-1, 1, 2}) {
          auto inds = vertex_inds(3, {dim_arr, sign_arr, r});
          for (int i_vert = 0; i_vert < 4; ++i_vert) {
            if (   &elem_arr[0]->shape().vertex(inds[0][i_vert])
                == &elem_arr[1]->shape().vertex(inds[1][i_vert])) {
              rotate = r;
            }
          }
        }
        _connect(elem_arr, {dim_arr, sign_arr, rotate});
      }
    }
    // rebuild `extrude_cons`
    for (int i_con = 0; i_con < (Int)def.cons.size(); ++i_con) {
      auto& con = def.cons[i_con];
      if (def.cons[i_con]) {
        if (def.cons[i_con]->element(1).tree && !def.cons[i_con]->element(0).tree) extrude_cons.push_back(con.get());
      }
    }
  }

  #pragma omp parallel for
  for (auto& vert : all_verts) {
    vert.record.clear();
  }
  purge();
  _offset_vertices(.01);
  #pragma omp parallel for
  for (Int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].active_shape().is_new = false;
  }
  {
    Task_message message(printers::info, "  Post-edge-matching mesh optimization", "\n", "  ");
    _optimize(1, 10, true);
  }
  Int n_failed = 0;
  for (auto& vert : all_verts) {
    n_failed += !vert.snap_to(_get_snapping_target(vert, vert.unwarped_point()));
  }
  if (n_failed) {
    printers::warn(format_str(200, "%li vertices could not be snapped to the surface.\n", n_failed), true);
  }
  #pragma omp parallel for
  for (auto& vert : all_verts) vert.record.clear();

  auto snap_intersections = [this](next::Boundary_block& block) {
    Array<double> interior {block.interior().reshaped({whatever, 3})};
    const Basis& b = block.basis();
    bool failed = false;
    int i_face = block.element()->boundary_face();
    Array<double> line_points({2, interior.shape()[0], 3});
    int block_nd = block.n_dim();
    int n_point = interior.shape()[0];
    for (int i_point = 0; i_point < n_point; ++i_point) {
      std::vector<int> coords(block_nd);
      for (int i_dim = 0; i_dim < block_nd; ++i_dim) {
        coords[i_dim] = math::row_coordinate(block_nd, b.row_size - 2, i_dim, i_point) + 1;
      }
      std::vector<int> elem_coords = block.element_coords(coords);
      line_points(1)(i_point).vector() = block.element()->point(elem_coords);
      elem_coords[i_face/2] = b.row_size - 1 - elem_coords[i_face/2];
      line_points(0)(i_point).vector() = block.element()->point(elem_coords);
    }
    for (int i_point = 0; i_point < n_point; ++i_point) {
      Mat<3> p0 = line_points(0)(i_point).vector();
      Mat<3> p1 = line_points(1)(i_point).vector();
      auto sects = surf_geom->intersections(resize(p0, params.n_dim), resize(p1, params.n_dim));
      double sect = huge;
      for (double s : sects) if (s > 0.) sect = std::min(sect, s);
      if (sect < 2.) {
        interior(i_point).vector() = p0*(1 - sect) + p1*sect;
      } else {
        failed = true;
      }
    }
    return failed;
  };

  auto check_elems = [this](next::Boundary_block& block) {
    bool failed = false;
    for (next::Element_shape* e : block.dependent_elements()) {
      next::Face* face = e->boundary_face_3d();
      if (face && face != &block) face->reset();
      Array<double> points = e->points();
      int nd = params.n_dim;
      int n_check_point = math::pow(_basis.row_size + 1, nd);
      Array<double> check_points({nd, n_check_point});
      Gauss_lobatto check_basis(_basis.row_size + 1);
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        check_points(i_dim).vector() = math::hypercube_matvec(_basis.interpolate(check_basis.nodes()),
                                                              points(i_dim).vector());
      }
      Array<double> jacobian({nd, nd, n_check_point});
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        for (int j_dim = 0; j_dim < nd; ++j_dim) {
          jacobian(i_dim)(j_dim).vector() = math::dimension_matvec(check_basis.diff_mat(),
                                                                   check_points(i_dim).vector(), j_dim);
        }
      }
      for (int i_check_point = 0; i_check_point < n_check_point; ++i_check_point) {
        Mat<dyn, dyn> point_jac(nd, nd);
        for (int i_dim = 0; i_dim < nd; ++i_dim) {
          for (int j_dim = 0; j_dim < nd; ++j_dim) {
            point_jac(i_dim, j_dim) = jacobian(i_dim)(j_dim)[i_check_point];
          }
        }
        if (!(point_jac.determinant() > 0)) failed = true;
      }
    }
    return failed;
  };

  auto plain_snap = [snap_intersections, check_elems](next::Boundary_block& block) {
    block.reset();
    bool failed = snap_intersections(block);
    if (!failed) failed = check_elems(block);
    if (failed) block.reset();
  };

  // snap edges to the surface (regardless of dimensionality)
  auto edges_2d = _blocks.edges_2d();
  #pragma omp parallel for
  for (auto& edge : edges_2d) plain_snap(edge);
  auto faces_3d = _blocks.faces_3d();
  #pragma omp parallel for
  for (auto& face : faces_3d) {
    for (int i_edge = 0; i_edge < 4; ++i_edge) face.edge(i_edge).reset();
  }
  #pragma omp parallel for
  for (auto& face : faces_3d) {
    for (int i_edge = 0; i_edge < 4; ++i_edge) {
      auto& edge = face.edge(i_edge);
      if (!edge.glued()) {
        std::vector<next::Element_shape*> dependent_elems = edge.dependent_elements();
        std::sort(dependent_elems.begin(), dependent_elems.end(), std::less());
        std::vector<std::unique_ptr<Lock::Set>> sets;
        for (next::Element_shape* e : dependent_elems) sets.emplace_back(new Lock::Set(e->lock));
        bool failed = false;
        if (edge.snapped_edge >= 0) {
          HEXED_ASSERT(edge.snapped_edge < edges.size(), "clearly erroneous `snapped_edge` value")
          auto& geom_edge = edges[edge.snapped_edge];
          Array<double> interior = edge.interior();
          for (int i_node = 0; i_node < interior.shape()[0]; ++i_node) {
            auto node = interior(i_node);
            auto nearest = geom_edge.nearest_point(node.vector());
            if (nearest.index >= 0) {
              node = geom_edge.nodes()(nearest.index);
            } else {
              failed = true;
            }
          }
        } else {
          failed = snap_intersections(edge);
        }
        if (!failed) failed = check_elems(edge);
        if (failed) edge.reset();
      }
    }
  }
  // snap face interiors (if 3D) to surface
  #pragma omp parallel for
  for (auto& face : faces_3d) plain_snap(face);
  ++_stopwatch["update"]["fit surface"].work_units_completed;
}

void Accessible_mesh::_optimize(int min_pow, int max_pow, bool check_snapping) {
  Stopwatch_tree::Starter sw_opt(_stopwatch["update"]["fit surface"]["optimization"]);
  auto verts = _blocks.verts();
  auto bverts = _blocks.boundary_verts();
  auto edges = surf_geom->edges();
  // determine which vertices are on extremal boundaries
  #pragma omp parallel for
  for (auto& vert : verts) {
    vert.record.resize(2*params.n_dim + 1);
    for (int i = 0; i < 2*params.n_dim + 1; ++i) vert.record[i] = 0;
  }
  #pragma omp parallel for
  for (auto& con : def.bound_cons) {
    int bc_sn = con->bound_cond_serial_n();
    if (bc_sn < 2*params.n_dim + 1) {
      std::vector<int> inds = vertex_inds(params.n_dim, con->get_direction())[0];
      for (int i_vert : inds) {
        #pragma omp atomic write
        con->element().active_shape().vertex(i_vert).record[bc_sn] = 1;
      }
    }
  }
  // boundary connections haven't been formed yet, so we have to do this manually
  #pragma omp parallel for
  for (auto& vert : bverts) {
    vert.record[2*params.n_dim] = 1;
  }
  History_monitor obj_monitor(.3, 100);
  History_monitor dist_monitor(.3, 100);
  double starting_objective = -1;
  double objective = 0;
  #if HEXED_VIS_MESH_OPT
  int id = rand()%1000;
  printers::info("id: " + std::to_string(id) + "\n");
  #endif
  std::vector<next::Vertex*> mobile_verts;
  for (auto& vert : verts) if (vert.mobile()) mobile_verts.push_back(&vert);
  auto compare = [](next::Vertex* vert0, next::Vertex* vert1) {return vert0->point({})[1] < vert1->point({})[1];};
  std::sort(mobile_verts.begin(), mobile_verts.end(), compare);
  #pragma omp parallel for
  for (auto vert : mobile_verts) vert->compute_depends();
  Int snaps_failed = 0;
  double total_dist = 0;
  Stopwatch watch;
  watch.start();
  double last_time = 0;
  std::string message;
  Int snap_succeeded = -1;
  next::Vertex::misses = 0;
  next::Vertex::tries = 0;
  for (Int i_relax = 0;
       i_relax < 1000 && (i_relax < 30
                          || (snaps_failed == 0 && obj_monitor.max() - obj_monitor.min()
                                                   > 1e-2*(std::abs(obj_monitor.max()) + std::abs(obj_monitor.min())))
                          || dist_monitor.max() - dist_monitor.min() > 1e-2*dist_monitor.min());
       ++i_relax) {
    #if HEXED_VIS_MESH_OPT
    {
      auto faces = _blocks.faces_3d();
      for (auto& f : faces) {
        for (int i_edge = 0; i_edge < 4; ++i_edge) f.edge(i_edge).reset();
      }
      for (auto& f : faces) {
        for (int i_edge = 0; i_edge < 4; ++i_edge) f.reset();
      }
      auto edges = _blocks.edges_2d();
      #pragma omp parallel for
      for (auto& e : edges) e.reset();
      visualize("default", "meshing_diagnostic" + std::to_string(id) + "_" + std::to_string(i_relax), (double)i_relax);
      std::string fname = "vertex_nearest" + to_string(id) + "_" + to_string(i_relax);
      auto vis = Visualizer::create("default", 3, 1, fname, {}, (double)i_relax, Visualizer::block);
      for (auto& vert : bverts) {
        Array<double> pos({3, 2});
        Mat<3> p0 = vert.unwarped_point();
        Mat<3> p1 = _get_snapping_target(vert, p0);
        for (int i_dim = 0; i_dim < 3; ++i_dim) {
          pos(i_dim)[0] = p0(i_dim);
          pos(i_dim)[1] = p1(i_dim);
        }
        vis->write_block(pos, Array<double>({}));
      }
    }
    #endif
    // snap vertices to surface boundary
    Mat<> o = tree->origin();
    double tns = tree->nominal_size();
    double objective_diff = 0;
    bool try_snap = true;
    if (try_snap) {
      snaps_failed = 0;
      total_dist = 0;
    }
    {
      Stopwatch_tree::Starter sw_relax(_stopwatch["update"]["fit surface"]["optimization"]["relaxation"]);
      Int n_chunk = mobile_verts.size()/config::n_threads + 1;
      #pragma omp parallel for reduction(+:objective_diff) reduction(+:snaps_failed) reduction(+:total_dist)
      for (int i_thread = 0; i_thread < config::n_threads; ++i_thread) {
      Int loop_max = std::min<Int>(mobile_verts.size(), (i_thread + 1)*n_chunk);
      for (int i_vert = i_thread*n_chunk; i_vert < loop_max; ++i_vert) {
        auto& vert = *mobile_verts[i_vert];
        auto satisfy = [&](Mat<3> p)->Mat<3> {
          for (int i_dim = 0; i_dim < (int)o.size(); ++i_dim) {
            p(i_dim) = std::max(p(i_dim), o(i_dim));
            p(i_dim) = std::min(p(i_dim), o(i_dim) + tns);
          }
          if (vert.record[2*params.n_dim]) {
            for (auto n : vert.neighbors()) if (n) {
              if ((int)n->record.size() == 2*params.n_dim + 1) {
                if (!n->record[2*params.n_dim]) {
                  auto seq = Eigen::seqN(0, params.n_dim);
                  Mat<> start = n->point({})(seq);
                  Mat<> end = p(Eigen::seqN(0, params.n_dim));
                  std::vector<double> intersections = surf_geom->intersections(start, end, false);
                  double min_sect = 1;
                  for (double s : intersections) min_sect = std::min(min_sect, s);
                  p(seq) = start + min_sect*(p(seq) - start);
                }
              }
            }
          }
          return p;
        };
        auto get_target = [&vert, this](Mat<3> p)->Mat<3>{return _get_snapping_target(vert, p);};
        bool on_surface = false;
        for (int i = 0; i < 2*params.n_dim + 1; ++i) on_surface = on_surface || vert.record[i];
        next::Vertex::Improve_quality_result iqr;
        if (on_surface) {
          iqr = vert.improve_quality(get_target, satisfy, true, try_snap);
        } else {
          iqr = vert.improve_quality();
        }
        objective_diff += iqr.objective_diff;
        snaps_failed += iqr.snap_failed;
        total_dist += iqr.target_dist;
      }
      }
      _stopwatch["update"]["fit surface"]["optimization"]["relaxation"].work_units_completed += mobile_verts.size();
    }
    double prev_obj = objective;
    objective = 0;
    {
      Stopwatch_tree::Starter sw_assess(_stopwatch["update"]["fit surface"]["optimization"]["assessment"]);
      #pragma omp parallel for reduction(+:objective)
      for (auto& vert : verts) objective += vert.quality_objective();
      _stopwatch["update"]["fit surface"]["optimization"]["assessment"].work_units_completed += verts.size();
    }
    if (starting_objective < 0) starting_objective = objective;
    if (i_relax) {
      if (std::abs(objective_diff - (objective - prev_obj)) > 1e-4*verts.size()) {
        printers::warn(" Warning: ", true);
        printers::warn(format_str(100, "inaccurate objective change: %e vs %e (please report as a bug)\n",
                                  -objective_diff, objective - prev_obj));
      }
    }
    obj_monitor.add_sample(i_relax, objective - starting_objective);
    dist_monitor.add_sample(i_relax, total_dist);
    message = format_str(
      400,
      "   "
      " Iteration = %4li;"
      " Objective = %.18e (%+.5e);"
      " Number of vertex snaps failed = %li;"
      " Total distance from surface = %.5e;"
      , i_relax, objective, objective - starting_objective, snaps_failed, total_dist
    );
    if (watch.time() > last_time) {
      last_time += .1;
      printers::info(message, false, true);
    }
    if (try_snap && snaps_failed == 0 && snap_succeeded < 0) snap_succeeded = i_relax;
    if (i_relax > 2*snap_succeeded && snap_succeeded >= 0 && i_relax > 30) break;
  }
  printers::info(message, false, true);
  printers::info("\n");
  printers::info(to_string(next::Vertex::misses) + " misses out of " + to_string(next::Vertex::tries) + " tries ("
                 + to_string(double(next::Vertex::misses)/next::Vertex::tries) + ")");
  ++_stopwatch["update"]["fit surface"]["optimization"].work_units_completed;
}

Storage_params incr_res_cache(Storage_params params) {
  params.n_stage += 1;
  return params;
}

Accessible_mesh::Accessible_mesh(Storage_params params_arg, double root_size_arg, Turbulence_model turb)
: params{params_arg}
, n_vert{math::pow(2, params.n_dim)}
, root_sz{root_size_arg}
, car{params, root_sz}
, def{incr_res_cache(params), root_sz}
, def_as_car{def.elements()}
, elems{car.elements(), def_as_car}
, kernel_elems{elems}
, elem_cons{car.element_connections(), def.element_connections()}
, bound_face_cons{car.bound_face_con_view, def.bound_face_con_view}
, bound_cons{car.boundary_connections(), def.boundary_connections()}
, def_face_cons{def.elem_face_con_v, bound_face_cons}
, ref_face_v{car.refined_faces(), def.refined_faces()}
, surf_bc_sn{-1} // set to -1 to prevent uninitialized comparisons
, verts_are_reset{false}
, _mask_levels{0}
, _basis(params.row_size)
, _blocks(params.n_dim, _basis)
, _n_verts{0}
, _stopwatch("mesh")
, _turb{turb}
, buffer_dist{std::sqrt(params.n_dim)/2} // if you're getting snapping problems, try multiplying this by 2
{
  def.face_con_v = def_face_cons;
  _stopwatch.emplace("update", "update");
  _stopwatch["update"].emplace("refinement", "refinement");
  _stopwatch["update"].emplace("extrusion", "extrusion");
  _stopwatch["update"].emplace("fit surface", "fit");
  _stopwatch["update"]["fit surface"].emplace("optimization", "optimization");
  _stopwatch["update"]["fit surface"]["optimization"].emplace("relaxation", "vertex update");
  _stopwatch["update"]["fit surface"]["optimization"].emplace("assessment", "vertex assessment");
  _stopwatch.work_units_completed = 1;
}

Accessible_mesh::~Accessible_mesh() {
  // delete connections before anything else so that deleting elements doesn't create dangling references
  car.purge_connections(criteria::always);
  def.purge_connections(criteria::always);
}

int Accessible_mesh::add_element(int ref_level, bool is_deformed, std::vector<Int> position, Mat<> origin, int aniso_ref_level, int surface_face) {
  int sn = container(is_deformed).emplace(ref_level, position, origin, aniso_ref_level);
  Element& elem = element(ref_level, is_deformed, sn);
  elem.create_shape(_blocks, surface_face);
  return sn;
}

int Accessible_mesh::add_element(int ref_level, bool is_deformed, std::vector<Int> position) {
  return add_element(ref_level, is_deformed, position, Mat<>::Zero(params.n_dim));
}

Element& Accessible_mesh::element(int ref_level, bool is_deformed, int serial_n) {
  return container(is_deformed).at(ref_level, serial_n);
}

void Accessible_mesh::_connect_shapes(Element& elem0, Element& elem1, Connection_direction dir) {
  next::Element_shape* shapes [2] {&elem0.shape(), &elem1.shape()};
  shapes[0]->connect(*shapes[1], dir);
  next::Element_shape* fake_shapes [2] {elem0.fake_shape(), elem1.fake_shape()};
  if (fake_shapes[0] || fake_shapes[1]) {
    for (int i = 0; i < 2; ++i) shapes[i] = fake_shapes[i] ? fake_shapes[i] : shapes[i];
    shapes[0]->connect(*shapes[1], dir);
  }
}

void Accessible_mesh::_connect(std::array<Element*, 2> el_ar, Con_dir<Element> direction) {
  car.cons.emplace_back(new Element_face_connection<Element>(el_ar, direction));
  _connect_shapes(*el_ar[0], *el_ar[1], direction);
}

void Accessible_mesh::_connect(std::array<Deformed_element*, 2> el_ar, Con_dir<Deformed_element> direction) {
  def.cons.emplace_back(new Element_face_connection<Deformed_element>(el_ar, direction));
  _connect_shapes(*el_ar[0], *el_ar[1], direction);
}

template <typename Elem_t>
void Accessible_mesh::_connect_shapes(Elem_t* coarse, std::vector<Elem_t*> fine, Con_dir<Deformed_element> dir,
              std::array<bool, 2> stretch) {
  std::vector<next::Element_shape*> fine_shapes;
  std::vector<next::Element_shape*> fake_fine_shapes;
  for (int i = 0; i < 1 + stretch[0]; ++i) {
    for (Elem_t* elem : fine) {
      for (int i = 0; i < 1 + stretch[1]; ++i) {
        fine_shapes.push_back(&elem->shape());
        fake_fine_shapes.push_back(elem->fake_shape());
      }
    }
  }
  coarse->shape().connect(fine_shapes, dir);
  HEXED_ASSERT(   std::all_of(fake_fine_shapes.begin(), fake_fine_shapes.end(), [](void* p)->bool{return  p;})
               || std::all_of(fake_fine_shapes.begin(), fake_fine_shapes.end(), [](void* p)->bool{return !p;}),
               "All of the fine elements must have fake shapes or none.");
  next::Element_shape* coarse_fake = coarse->fake_shape();
  if (coarse_fake || fake_fine_shapes[0]) {
    next::Element_shape* shape = coarse_fake ? coarse_fake : &coarse->shape();
    shape->connect(fake_fine_shapes[0] ? fake_fine_shapes : fine_shapes, dir);
  }
}

void Accessible_mesh::_connect(Element* coarse, std::vector<Element*> fine, Con_dir<Deformed_element> dir) {
  HEXED_ASSERT(dir.i_dim[0] == dir.i_dim[1], "dimensions in Cartesian hanging-node connection must match");
  car.ref_face_cons[params.n_dim - 1].emplace_back(
    new Refined_connection<Element>(coarse, fine, {dir.i_dim[0]}, dir.face_sign[1])
  );
  _connect_shapes(coarse, fine, dir, {false, false});
}

void Accessible_mesh::_connect(Deformed_element* coarse, std::vector<Deformed_element*> fine,
                               Con_dir<Deformed_element> dir, std::array<bool, 2> stretch) {
  def.ref_face_cons[math::log(2, fine.size())].emplace_back(
    new Refined_connection<Deformed_element>(coarse, fine, dir, false, stretch)
  );
  _connect_shapes(coarse, fine, dir, stretch);
}

void Accessible_mesh::connect_cartesian(int ref_level, std::array<Int, 2> serial_n, Con_dir<Element> direction,
                                        std::array<bool, 2> is_deformed) {
  std::array<Element*, 2> el_ar;
  for (int i_side : {0, 1}) el_ar[i_side] = &element(ref_level, is_deformed[i_side], serial_n[i_side]);
  _connect(el_ar, direction);
}

void Accessible_mesh::connect_deformed(int ref_level, std::array<Int, 2> serial_n,
                                       Con_dir<Deformed_element> direction) {
  if ((direction.i_dim[0] == direction.i_dim[1]) && (direction.face_sign[0] == direction.face_sign[1])) {
    throw std::runtime_error("attempt to connect faces of same sign along same dimension which is forbidden");
  }
  std::array<Deformed_element*, 2> el_ar;
  for (int i_side : {0, 1}) {
    el_ar[i_side] = &def.elems.at(ref_level, serial_n[i_side]);
  }
  _connect(el_ar, direction);
}

void Accessible_mesh::connect_hanging(int coarse_ref_level, Int coarse_serial, std::vector<Int> fine_serial,
                                      Con_dir<Deformed_element> dir, bool coarse_deformed,
                                      std::vector<bool> fine_deformed, std::array<bool, 2> stretch) {
  bool is_car = !coarse_deformed;
  for (bool fine_def : fine_deformed) is_car = (is_car||!fine_def);
  if (is_car) {
    Element* coarse = &element(coarse_ref_level, coarse_deformed, coarse_serial);
    std::vector<Element*> fine;
    for (int i_fine = 0; i_fine < n_vert/2; ++i_fine) {
      fine.push_back(&element(coarse_ref_level + 1, fine_deformed[i_fine], fine_serial[i_fine]));
    }
    HEXED_ASSERT((dir.i_dim[0] == dir.i_dim[1]) && (dir.face_sign[0] != dir.face_sign[1]),
                 "attempted to form a cartesian hanging-node connection with incompatible `Con_dir`.");
    _connect(coarse, fine, dir);
  } else {
    Deformed_element* coarse = &def.elems.at(coarse_ref_level, coarse_serial);
    std::vector<Deformed_element*> fine;
    for (unsigned i_fine = 0; i_fine < fine_serial.size(); ++i_fine) {
      fine.push_back(&def.elems.at(coarse_ref_level + 1, fine_serial[i_fine]));
    }
    _connect(coarse, fine, dir, stretch);
  }
}

int Accessible_mesh::add_boundary_condition(Flow_bc* flow_bc) {
  bound_conds.emplace_back(flow_bc);
  // no reason to delete boundary conditions, so the serial number can just be the index
  return bound_conds.size() - 1;
}

void Accessible_mesh::connect_boundary(int ref_level, bool is_deformed, Int element_serial_n, int i_dim, int face_sign, int bc_serial_n) {
  // create boundary condition
  HEXED_ASSERT(bc_serial_n < int(bound_conds.size()), "demand for non-existent `Boundary_condition`");
  Flow_bc& bc = *bound_conds[bc_serial_n];
  if (is_deformed) {
    def.bound_cons.emplace_back(new Typed_bound_connection<Deformed_element>(
      def.elems.at(ref_level, element_serial_n), i_dim, face_sign, bc_serial_n, bc.n_prescribed(params.n_dim)
     ));
  } else {
    car.bound_cons.emplace_back(new Typed_bound_connection<Element>(
      car.elems.at(ref_level, element_serial_n), i_dim, face_sign, bc_serial_n, bc.n_prescribed(params.n_dim)
    ));
  }
}

void Accessible_mesh::disconnect_boundary(int bc_sn) {
  erase_if(car.bound_cons, [bc_sn](std::unique_ptr<Typed_bound_connection<         Element>>& con){return con->bound_cond_serial_n() == bc_sn;});
  erase_if(def.bound_cons, [bc_sn](std::unique_ptr<Typed_bound_connection<Deformed_element>>& con){return con->bound_cond_serial_n() == bc_sn;});
}

void Accessible_mesh::cleanup() {
  purge();
}

Mesh::Connection_validity Accessible_mesh::valid() {
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

void Accessible_mesh::assert_valid() {
  auto v = valid();
  VIS_ASSERT(v.n_redundant == 0 && v.n_missing == 0, format_str(200,
    "invalid mesh connectivity with %i redundancies and %i unconnected faces",
    v.n_redundant, v.n_missing
  ));
}

//! \cond helper classes and functions for Accessible_mesh::extrude
struct Empty_face {
  Deformed_element& elem;
  int i_dim;
  int face_sign;
};
struct Connection_plan {
  int ref_level;
  std::array<Int, 2> serial_ns;
  Con_dir<Deformed_element> dir;
};
struct Refined_connection_plan {
  int coarse_ref;
  Int coarse_sn;
  std::vector<Int> fine_sn;
  Con_dir<Deformed_element> dir;
  std::array<bool, 2> stretch;
};
bool aligned_same_dim(Con_dir<Deformed_element> dir, std::array<Int, 2> extrude_dim) {
  return (dir.i_dim[0] == dir.i_dim[1]) && (dir.face_sign[0] != dir.face_sign[1])
         && (extrude_dim[0] == extrude_dim[1]);
}
bool aligned_different_dim(Con_dir<Deformed_element> dir, std::array<Int, 2> extrude_dim) {
  return (dir.i_dim[0] == extrude_dim[1]) && (dir.i_dim[1] == extrude_dim[0]);
}
//! \endcond

void request_connection(Element& elem, int n_dim, int i_dim, bool i_sign, int j_dim, bool j_sign) {
  // record data at vertex which is on the face to be connected, on the face which was extruded from,
  // and if applicable has the minimum index to satisfy the above consitions.
  int i_vert =   j_sign*math::pow(2, n_dim - 1 - j_dim)
               + (1 - i_sign)*math::pow(2, n_dim - 1 - i_dim);
  auto& record = elem.shape().vertex(i_vert).record;
  // which element it is
  record.push_back(elem.refinement_level());
  record.push_back(elem.record); // `elem.record` = serial number
  // which face needs to be connected
  record.push_back(2*j_dim + j_sign);
  // extrusion direction for deciding which face connections are valid
  record.push_back(2*i_dim + i_sign);
}

void Accessible_mesh::extrude(bool collapse, double offset, bool force) {
  Stopwatch_tree::Starter sw_extrude(_stopwatch["update"]["extrusion"]);
  const int nd = params.n_dim;
  const int n_faces = 2*nd;
  { // initialize number of connections of each face to 0
    auto& elems = elements();
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      for (int i_face = 0; i_face < n_faces; ++i_face) {
        elems[i_elem].face_record[i_face] = 0;
      }
    }
  }
  { // initialize vertex records to empty
    auto verts = _blocks.verts();
    for (Int i_vert = 0; i_vert < verts.size(); ++i_vert) {
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

  // request connections for existing extruded elements
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
    auto verts = _blocks.verts();
    for (Int i_vert = 0; i_vert < verts.size(); ++i_vert) {
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
  for (auto face : empty_faces) {
    auto nom_pos = face.elem.nominal_position();
    nom_pos[face.i_dim] += 2*face.face_sign - 1;
    const int ref_level = face.elem.refinement_level();
    int sn = add_element(ref_level, true, nom_pos, face.elem.origin, face.elem.aniso_ref_level() + 1, 2*face.i_dim + face.face_sign);
    Con_dir<Deformed_element> dir {{face.i_dim, face.i_dim}, {!face.face_sign, bool(face.face_sign)}};
    auto& elem = def.elems.at(ref_level, sn);
    if (face.elem.fake_shape()) elem.split_shape(_blocks, face.elem, offset, 2*face.i_dim + face.face_sign);
    else elem.create_fake(_blocks);
    elem.record = sn;
    elem.needs_snapping = !force;
    if (collapse) {
      int stride = math::pow(2, nd - 1 - face.i_dim);
      for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
        int i_collapse = i_vert + (face.face_sign - (i_vert/stride)%2)*stride;
        elem.fake_shape()->vertex(i_vert).point({}) = face.elem.shape().vertex(i_collapse).point({});
      }
    }
    elem.fake_shape()->extruded_direction = 2*face.i_dim + face.face_sign;
    std::array<Deformed_element*, 2> el_arr {&elem, &face.elem};
    _connect(el_arr, dir);
    extrude_cons.push_back(def.cons.back().get());
    // record the faces that still need to be connected at a vertex which is guaranteed to be shared with prospective neighbors
    for (int j_dim = face.i_dim + 1; j_dim%nd != face.i_dim; ++j_dim) {
      j_dim = j_dim%nd;
      for (int face_sign = 0; face_sign < 2; ++face_sign) {
        int face_rec = face.elem.face_record[2*j_dim + face_sign];
        if (face_rec >= 2) {
          // if parent element has boundary connections on other faces
          def.bound_cons.emplace_back(new Typed_bound_connection<Deformed_element>(
            elem, j_dim, face_sign, face_rec - 2, bound_conds[face_rec - 2]->n_prescribed(nd)
          ));
        } else request_connection(elem, nd, face.i_dim, face.face_sign, j_dim, face_sign);
      }
    }
    if (offset > 0) {
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
    auto verts = _blocks.verts();
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
      for (bool (*aligned)(Con_dir<Deformed_element>, std::array<Int, 2>) : {&aligned_same_dim, &aligned_different_dim}) { \
        /* iterate through every possible pair of records created by an extruded elements above */ \
        for (int i_record = 0; i_record < int(vert.record.size()); i_record += n_record) { \
          for (int j_record = i_record + n_record; j_record < int(vert.record.size()); j_record += n_record) { \
            int ref_level = vert.record[i_record]; \
            Con_dir<Deformed_element> dir({ int(vert.record[i_record + 2]/2),  int(vert.record[j_record + 2]/2)}, \
                                          {bool(vert.record[i_record + 2]%2), bool(vert.record[j_record + 2]%2)}); \
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
      std::array<Int, 2> sn {vert.record[i_record + 1], vert.record[j_record + 1]}; \
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
      std::vector<Int> sn(n_vert/4); \
      Int* vert_rec = vert.record.data() + rec[which_fine]; \
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
        next::Vertex& fine_vert = elem.shape().vertex(iv); \
        /* vertex should have exactly one record */ \
        HEXED_ASSERT(fine_vert.record.size() == n_record, \
                     format_str(1000, "`fine_vert.record.size() == %li != n_record == %li` (position = [%e, %e, %e])", \
                                fine_vert.record.size(), n_record, \
                                fine_vert.point({})[0], fine_vert.point({})[1], fine_vert.point({})[2]).c_str()); \
        sn[1] = fine_vert.record[1]; \
        stretch_dim = extr_dim > free_dim; \
        fine_vert.record.erase(fine_vert.record.begin(), fine_vert.record.begin() + n_record); \
      } \
      std::array<bool, 2> stretch {false, false}; \
      stretch[stretch_dim] = true; \
      Int* coarse_rec = vert.record.data() + rec[!which_fine]; \
      ref_con_plans.push_back({(int)coarse_rec[0], coarse_rec[1], sn, {{(int)coarse_rec[2]/2, face_dim}, {bool(coarse_rec[2]%2), bool(face_sign)}}, stretch}); \
      /* erase unconnected face record to prevent duplicate connections */ \
      vert.record.erase(vert.record.begin() + j_record, vert.record.begin() + j_record + n_record); \
      vert.record.erase(vert.record.begin() + i_record, vert.record.begin() + i_record + n_record); \
      i_record -= n_record; /* move index to account for erased elements */ \
      break; /* since we found a match, we can move on to the next `i_record` */ \
    } \

  {
    auto verts = _blocks.verts();
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
    auto verts = _blocks.verts(); // note: need to rebuild vertex vector because face connections above have `eat`en vertices
    VERTEX_LOOP(CONNECT_SAME_CONFORMING)
    for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
      HEXED_ASSERT(verts[i_vert].record.empty(), format_str(100, "vertex detected with %lu unprocessed connection requests",
                                                            verts[i_vert].record.size()/n_record).c_str());
    }
  }
  for (auto con_plan : con_plans) {
    connect_deformed(con_plan.ref_level, con_plan.serial_ns, con_plan.dir);
  }
  _n_verts = _blocks.verts().size();
  ++_stopwatch["update"]["extrusion"].work_units_completed;
}

void Accessible_mesh::connect_rest(int bc_sn) {
  HEXED_ASSERT((Int)bound_conds.size() > bc_sn, "nonexistant boundary condition");
  HEXED_ASSERT(bound_conds[bc_sn], "BC pointer is null");
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
  car.connect_empty(bc_sn, *bound_conds[bc_sn]);
  def.connect_empty(bc_sn, *bound_conds[bc_sn]);
}

std::vector<Mesh::elem_handle> Accessible_mesh::elem_handles() {
  std::vector<Mesh::elem_handle> handles;
  for (bool is_deformed : {0, 1}) {
    for (auto ref_sn : container(is_deformed).elem_handles()) {
      handles.push_back({ref_sn[0], is_deformed, ref_sn[1]});
    }
  }
  return handles;
}

Element& Accessible_mesh::add_elem(bool is_deformed, Tree& t) {
  auto np = t.coordinates();
  int sn = add_element(t.refinement_level(), is_deformed, std::vector<Int>(np.begin(), np.end()), t.origin());
  auto& elem = element(t.refinement_level(), is_deformed, sn);
  elem.record = sn; // put the serial number in the record so it can be used for connections
  elem.tree.pair(t.elem);
  if (is_deformed) {
    t.def_elem = &def.elems.at(t.refinement_level(), sn);
  }
  return elem;
}

void Accessible_mesh::create_tree(std::vector<Flow_bc*> extremal_bcs, Mat<> origin) {
  // take ownership of bcs (do this first to avoid memory leak)
  std::vector<int> new_tree_bcs;
  //! \todo this could, in theory, be a resource leak because these are never erased if an exception is thrown...
  for (Flow_bc* fbc : extremal_bcs) new_tree_bcs.push_back(add_boundary_condition(fbc));
  HEXED_ASSERT(int(extremal_bcs.size()) == 2*params.n_dim, "`extremal_bcs` has wrong number of elements");
  HEXED_ASSERT(!tree, "each `Mesh` may only contain one tree");
  // add the tree
  tree_bcs = new_tree_bcs;
  tree.reset(new Tree(params.n_dim, root_sz, origin));
}

void Accessible_mesh::add_tree(std::vector<Flow_bc*> extremal_bcs, Mat<> origin) {
  create_tree(extremal_bcs, origin);
  auto& elem = add_elem(false, *tree);
  int sn = elem.record;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      connect_boundary(0, 0, sn, i_dim, sign, tree_bcs[2*i_dim + sign]);
    }
  }
  _n_verts = _blocks.verts().size();
}

bool Accessible_mesh::intersects_surface(Tree* t) {
  if (!surf_geom) return false;
  Mat<> center = t->nominal_position() + t->nominal_size()/2*Mat<>::Ones(params.n_dim);
  return !surf_geom->nearest_point(center, buffer_dist*t->nominal_size()).empty();
}

bool Accessible_mesh::is_surface(Tree* t) {
  return t->get_status() == 0;
}

void Accessible_mesh::set_surface(Surface_geom* geometry, Flow_bc* surface_bc, Eigen::VectorXd flood_fill_start) {
  Task_message tm0(printers::info, "  incorporating surface geometry", "\n");
  // take ownership of the surface geometries (do this first to avoid memory leak)
  surf_bc_sn = add_boundary_condition(surface_bc);
  surf_geom.reset(geometry);
  Int n_edges = surf_geom->edges().size();
  matched_vertices.clear();
  matched_vertices.resize(n_edges);
  matched_edges.clear();
  matched_edges.resize(n_edges);
  point_matched_vertices.clear();
  point_matched_vertices.resize(surf_geom->points().size());
  if (!tree) return;
  // identify surface elements
  auto& elems = elements();
  {
    Task_message tm1(printers::info, "    determining inside/outside", "\n");
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      if (elem.tree) if (intersects_surface(elem.tree.get())) elem.tree->set_status(0);
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
  }
  delete_bad_extrusions();
  deform();
  purge();
  connect_new<         Element>(0);
  connect_new<Deformed_element>(0);
  extrude(true);
  connect_rest(surf_bc_sn);
  _fit_surface();
  connect_rest(surf_bc_sn);
  _n_verts = _blocks.verts().size();
}

void Accessible_mesh::set_unref_locks(std::function<bool(Element&)> lock_if) {
  auto& elems = elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].unrefinement_locked = lock_if(elems[i_elem]);
  }
}

// forms connections for new tree elements of type `element_t` starting with the `starting_at`th one
template<typename element_t>
void Accessible_mesh::connect_new(int start_at) {
  auto& m = mbt<element_t>();
  auto elems = m.elems.elements();
  int nd = params.n_dim;
  // helper function for connecting refined elements
  auto connect_refined = [&](Element& elem, int i_dim, int sign, std::vector<Tree*> neighbors) {
    VIS_ASSERT(int(neighbors.size()) == math::pow(2, nd - 1),
               format_str(100, "bad number of neighbors %lu (thanks for nothing, ref level smoother)",
                          neighbors.size()))
    bool is_def = elem.get_is_deformed();
    for (Tree* neighbor : neighbors) {
      VIS_ASSERT(neighbor->elem, "hanging-node connection with nonexistant elements");
      is_def = is_def && neighbor->elem->get_is_deformed();
    }
    Con_dir<Deformed_element> dir {{i_dim, i_dim}, {!sign, bool(sign)}};
    if (is_def) {
      std::vector<Deformed_element*> fine;
      for (Tree* neighbor : neighbors) fine.push_back(neighbor->def_elem);
      _connect(elem.tree->def_elem, fine, dir);
    } else {
      std::vector<Element*> fine;
      for (Tree* neighbor : neighbors) fine.push_back(neighbor->elem.get());
      _connect(&elem, fine, dir);
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
            // if this element is at the boundary of the tree (as opposed to a surface geometry boundary),
            // set an extremal boundary condition
            if (neighbors.empty()) {
              m.bound_cons.emplace_back(new Typed_bound_connection<element_t>(
                elem, i_dim, sign, tree_bcs[2*i_dim + sign], bound_conds[tree_bcs[2*i_dim + sign]]->n_prescribed(nd)
              ));
            }
            // otherwise, if the element has not only a tree neighbor but also an element neighbor...
            else if (neighbors[0]->elem) {
              if (neighbors.size() == 1) {
                auto& other = *neighbors[0]->elem;
                // if ref levels are the same, make a conformal connection
                if (other.refinement_level() == elem.refinement_level()) {
                  if (elem.deformed() && other.deformed()) {
                    std::array<Deformed_element*, 2> el_ar {elem.tree->def_elem, neighbors[0]->def_elem};
                    _connect(el_ar, Con_dir<Deformed_element>{{i_dim, i_dim}, {bool(sign), !sign}});
                  } else {
                    std::array<Element*, 2> el_ar;
                    el_ar[!sign] = &elem;
                    el_ar[sign] = &other;
                    _connect(el_ar, Con_dir<Element>{i_dim});
                  }
                } else {
                  // if neighbor is coarser, form a hanging node connection
                  // but only if this is the fine element with the lowest coordinates,
                  // to prevent redundant connections from all the fine elements
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

void Accessible_mesh::refine_set_status(Tree* t) {
  t->refine();
  for (Tree* child : t->children()) {
    child->set_status(intersects_surface(child) - 1);
  }
}

// performs the actual refinement for all elements where the record has been set to 1
void Accessible_mesh::refine_by_record(bool is_deformed, int start, int end) {
  auto& elems = container(is_deformed).element_view();
  for (int i_elem = start; i_elem < end; ++i_elem) {
    auto& elem = elems[i_elem];
    if (elem.tree) {
      if (elem.record == 1) {
        elem.record = 2;
        refine_set_status(elem.tree.get());
        for (Tree* child : elem.tree->children()) {
          add_elem(is_deformed, *child).record = 0;
        }
      }
    }
  }
}

bool exists(Tree* tree) {
  if (tree) if (tree->elem) {
    int record;
    #pragma omp atomic read
    record = tree->elem->record;
    if (record != 2) return true;
  }
  return false;
}

// does this tree element need to be refined to satisfy ref level smoothness
bool Accessible_mesh::needs_refine(Tree* t) {
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

bool has_existent_children(Tree* t) {
  bool has = t->elem;
  for (Tree* child : t->children()) has = has || has_existent_children(child);
  return has;
}

// whether an element is currently set to be deformed at the end of the `update` sweep
bool is_def(Element& elem) {
  return (elem.get_is_deformed() && elem.record == 0) || (!elem.get_is_deformed() && elem.record == 3);
}

// delete elements that would create pathological extrusion topology
void Accessible_mesh::delete_bad_extrusions() {
  int nd = params.n_dim;
  auto& elems = elements();
  bool changed;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    if (elems[i_elem].has_shape()) elems[i_elem].active_shape().record = 0;
  }
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
    // delete elements that would create pathological offset geometry
    auto verts = _blocks.verts();
    for (auto& vert : verts) vert.record.clear();
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      if (!elem.tree || !elem.has_shape() || elem.record == 2) continue;
      for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
        for (bool sign : {0, 1}) {
          if (!exists(elem.tree->find_neighbor(math::sign(sign)*Eigen::VectorXi::Unit(params.n_dim, i_dim)))) {
            auto inds = vertex_inds(params.n_dim, {{i_dim, i_dim}, {sign, !sign}})[0];
            for (int i_vert : inds) {
              elem.active_shape().vertex(i_vert).record.push_back(2*i_dim + sign);
            }
          }
        }
      }
    }
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      if (!elem.tree || !elem.has_shape() || elem.record == 2) continue;
      for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
        auto& vert = elem.active_shape().vertex(i_vert);
        for (int i = 0; i < Int(vert.record.size()); ++i) {
          for (int j = 0; j < Int(vert.record.size()); ++j) if (i != j) {
            if (vert.record[i]/2 == vert.record[j]/2 && vert.record[i]%2 != vert.record[j]%2) {
              changed = true;
              elem.record = 2;
              vert.record.clear();
            }
          }
        }
      }
    }
  } while (changed);
  auto verts = _blocks.verts();
  for (auto& vert : verts) vert.record.clear();
}

void Accessible_mesh::deform() {
  auto& elems = elements();
  int nd = params.n_dim;
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

void Accessible_mesh::purge() {
  if (tree) {
    // delete obsolete elements of `extrude_cons`
    erase_if(extrude_cons, [](Element_face_connection<Deformed_element>* con){return con->element(0).record == 2 || con->element(1).record == 2;});
    // delete connections to old elements (has to happen before deleting elements or else use after free)
    car.purge_connections();
    def.purge_connections();
    // delete old elements
    car.elems.purge();
    def.elems.purge();
    // delete old matched vertices and edges
    _blocks.verts(); // evaluating `verts` and `boundary_sides` automatically purges the vertex and face/edge lists
    _blocks.boundary_sides();
    std::erase(point_matched_vertices, nullptr);
    for (auto& vec : matched_vertices) std::erase(vec, nullptr);
    for (auto& vec : matched_edges) std::erase(vec, nullptr);
  }
}

bool Accessible_mesh::update(std::function<bool(Element&)> refine_criterion,
                             std::function<bool(Element&)> unrefine_criterion) {
  /* `Element::record` is used to identify which elements are going to be modified.
   * 0 => do nothing
   * 1 => refine
   * -1 => unrefine
   * 2 => delete
   * 3 => toggle deformity
   */
  HEXED_ASSERT(tree, "need a tree to refine");
  Stopwatch_tree::Starter sw_update(_stopwatch["update"]);
  int nd = params.n_dim;
  auto& elems = elements();
  Int n_before = elems.size();
  {
    Stopwatch_tree::Starter sw_refine(_stopwatch["update"]["refinement"]);
    // decide which elements to (un)refine
    #pragma omp parallel for // parallelize this part since `predicate` could be expensive
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      elem.record = 0;
      bool ref = refine_criterion(elem);
      bool unref = unrefine_criterion(elem);
      if (ref && !unref) elem.record = 1;
      else if (unref && !ref) elem.record = -1;
    }
    // pass refinement requests of extruded elements to their extrusion parents
    #pragma omp parallel for
    for (auto con : extrude_cons) {
      auto& inside = con->element(1);
      Lock::Set s(inside.lock);
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
        if (elem.tree && elem.record != 2) if (is_surface(elem.tree.get())) elem.record = 2;
      }
    }
    // incremental flood fill
    do {
      changed = false;
      // synchronize refinement level of surface elements with their non-surface neighbors
      for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
        auto& elem = elems[i_elem];
        if (exists(elem.tree.get())) {
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
          if (exists(elem.tree.get())) {
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
          if (needs_refine(elem.tree.get())) {
            changed = true;
            elem.record = 1;
          }
        }
      }
      for (bool is_deformed : {0, 1}) refine_by_record(is_deformed, 0, container(is_deformed).element_view().size());
    } while (changed);
    // set extruded elements to be deleted
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      if (!elems[i_elem].tree) elems[i_elem].record = 2;
    }
    purge();
    delete_bad_extrusions();
    purge();
    // connect new elements
    connect_new<         Element>(0);
    connect_new<Deformed_element>(0);
    deform();
    purge();
    connect_new<         Element>(0);
    connect_new<Deformed_element>(0);
    ++_stopwatch["update"]["refinement"].work_units_completed;
  }
  extrude(true);
  if (surf_geom) {
    connect_rest(surf_bc_sn);
    for (auto& ptr : point_matched_vertices) ptr.set();
    for (Int i_edge = 0; i_edge < (Int)surf_geom->edges().size(); ++i_edge) {
      matched_vertices[i_edge].clear();
      matched_edges[i_edge].clear();
    }
    _fit_surface();
    connect_rest(surf_bc_sn);
    _n_verts = _blocks.verts().size();
  }
  _n_verts = _blocks.verts().size();
  _stopwatch["update"].work_units_completed += 1;
  Int n_after = elems.size();
  printers::info(format_str(100, "  %li net elements created\n", n_after - n_before));
  for (auto& con : car.cons) {
    HEXED_ASSERT(!con->element(0).deformed() || !con->element(1).deformed(), "cartesian connection between deformed elements");
  }
  return n_after - n_before;
}

void update_pos(next::Vertex& vert, Mat<3> pos) {
  if ((pos - vert.point({})).norm() < vert.nominal_size()) vert.set_pos(pos);
}

Accessible_mesh::Masked_mesh::Masked_mesh(Accessible_mesh& mesh, const Basis& basis,
                                          std::function<bool(Element&)> mask)
: kernel_mesh {
    mesh.params.n_dim,
    mesh.params.row_size,
    mesh.params.n_var,
    mesh._mask_levels,
    basis,
    mesh._turb,
    _masked_car_cons.slice,
    _masked_def_cons.slice,
    _masked_car_elems.slice,
    _masked_def_elems.slice,
    _masked_elems.slice,
    _masked_ref_faces.slice,
  },
  bound_cons{_masked_bound_cons.slice}
{
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < mesh.elems.size(); ++i_elem) {
    if (mesh.elems[i_elem]._mask >= mesh._mask_levels - 1 && mask(mesh.elems[i_elem])) {
      mesh.elems[i_elem]._mask = mesh._mask_levels;
    }
  }
  #define MASK_REF_CONS(mbt) \
    for (int i_con = 0; i_con < mbt.refined_connections().size(); ++i_con) { \
      auto& con = mbt.refined_connections()[i_con]; \
      con.refined_face.coarse_mask = con.coarse_element().mask(); \
      for (int i_fine = 0; i_fine < con.n_fine_elements(); ++i_fine) { \
        con.refined_face.fine_masks[i_fine] = con.connection(i_fine).element(!con.order_reversed()).mask(); \
      } \
      for (int i_fine = con.n_fine_elements(); i_fine < 4; ++i_fine) con.refined_face.fine_masks[i_fine] = 0; \
    }
  #pragma omp parallel for
  MASK_REF_CONS(mesh.car)
  #pragma omp parallel for
  MASK_REF_CONS(mesh.def)
  #undef MASK_REF_CONS
  _masked_elems.populate(mesh.elems, [&](Element& elem){return elem.mask() >= mesh._mask_levels;});
  _masked_car_elems.populate(mesh.car.elements(), [&](Element& elem){return elem.mask() >= mesh._mask_levels;});
  _masked_def_elems.populate(mesh.def.elements(), [&](Element& elem){return elem.mask() >= mesh._mask_levels;});
  _masked_car_cons.populate(mesh.car.kernel_connections(), [&](Kernel_connection& con){return con.mask() >= mesh._mask_levels;});
  _masked_def_cons.populate(mesh.def.kernel_connections(), [&](Kernel_connection& con){return con.mask() >= mesh._mask_levels;});
  _masked_ref_faces.populate(mesh.ref_face_v, [&](Refined_face& face){return face.mask() >= mesh._mask_levels;});
  _masked_bound_cons.populate(mesh.bound_cons, [&](Boundary_connection& con){return con.mask() >= mesh._mask_levels;});
  ++mesh._mask_levels;
}

void Accessible_mesh::reset_masks() {
  _mask_levels = 0;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem]._mask = 0;
  }
}

std::vector<std::unique_ptr<Accessible_mesh::Masked_mesh>> Accessible_mesh::preti_masks(const Basis& basis) {
  reset_masks();
  std::vector<std::unique_ptr<Masked_mesh>> masks;
  masks.emplace_back(new Masked_mesh(*this, basis));
  while (masks.back()->kernel_mesh.elems.size()) {
    masks.emplace_back(new Masked_mesh(*this, basis, [this](Element& elem){return elem.aniso_ref_level() >= _mask_levels;}));
  }
  masks.pop_back();
  return masks;
}

void Accessible_mesh::reset_verts() {
  #if 0
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
        Lock::Set s(vert.lock);
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
      Lock::Set s(vert.lock);
      int face_sign = (i_vert/stride)%2;
      if (face_sign == dir.face_sign[1]) {
        vert.pos = elem.vertex(i_vert + stride*(dir.face_sign[0] - face_sign)).pos;
      }
    }
  }
  verts_are_reset = true;
  #endif
}

void Accessible_mesh::restore_verts() {
  #if 0
  verts_are_reset = false;
  auto verts = vertices();
  #pragma omp parallel for
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
    verts[i_vert].pos = verts[i_vert].temp_vector;
  }
  #endif
}

template <typename T>
void h5_write_row(H5::DataSet& dset, int cols, int i_row, T* data) {
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
void h5_write_value(H5::DataSet& dset, int i_row, T data) {
  h5_write_row(dset, 1, i_row, &data);
}

template <typename T>
void h5_add_attr(H5::H5Object& obj, std::string name, T value, H5::DataType dtype = H5::PredType::NATIVE_INT) {
  hsize_t attr_dim = 1;
  H5::DataSpace dspace(1, &attr_dim);
  auto attr = obj.createAttribute(name, dtype, dspace);
  attr.write(dtype, &value);
}

template <typename T = int>
T h5_get_attr(H5::H5Object& obj, std::string name, H5::DataType dtype = H5::PredType::NATIVE_INT) {
  auto attr = obj.openAttribute(name.c_str());
  T value;
  attr.read(dtype, &value);
  return value;
}

Storage_params read_params(std::string file_name) {
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

double read_root_sz(std::string file_name) {
  H5::H5File file(file_name + ".mesh.h5", H5F_ACC_RDONLY);
  return h5_get_attr<double>(file, "root_size", H5::PredType::NATIVE_DOUBLE);
}

template <typename T>
void h5_read_row(H5::DataSet& dset, int cols, int i_row, T* data) {
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
T h5_read_value(H5::DataSet& dset, int i_row) {
  T data;
  h5_read_row(dset, 1, i_row, &data);
  return data;
}

void Accessible_mesh::write(std::string name) {
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
  auto verts = _blocks.verts();
  dims[0] = verts.size();
  dims[1] = 3;
  auto dset = file.createDataSet("vertices/position", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, dims));
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
    auto& vert = verts[i_vert];
    Mat<3> p = vert.point({});
    h5_write_row(dset, 3, i_vert, p.data());
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
  auto aniso_ref_level_dset = file.createDataSet("/elements/aniso_refinement_level", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    int vert_inds [8] {};
    for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
      vert_inds[i_vert] = elem.shape().vertex(i_vert).record[0];
    }
    h5_write_row(vert_dset, params.n_vertices(), i_elem, vert_inds);
    std::vector<int> nom_pos;
    for (Int p : elem.nominal_position()) nom_pos.push_back(p);
    h5_write_row(nom_pos_dset, params.n_dim, i_elem, nom_pos.data());
    h5_write_value(is_def_dset, i_elem, elem.get_is_deformed());
    h5_write_value(ref_level_dset, i_elem, elem.refinement_level());
    h5_write_value(aniso_ref_level_dset, i_elem, elem.aniso_ref_level());
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
  dims[1] = 11;
  auto ref_con_dset = file.createDataSet("/connections/refined", H5::PredType::NATIVE_INT, H5::DataSpace(2, dims));
  #define WRITE_REF_CONS(start, cons) \
    for (int i_con = 0; i_con < cons.size(); ++i_con) { \
      auto& con = cons[i_con]; \
      int data[11] {}; \
      data[0] = con.coarse_element().record; \
      for (int i_fine = 0; i_fine < con.n_fine_elements(); ++i_fine) { \
        data[1 + i_fine] = con.connection(i_fine).element(!con.order_reversed()).record; \
      } \
      for (int i_fine = con.n_fine_elements(); i_fine < 4; ++i_fine) data[1 + i_fine] = -1; \
      for (int i_dim : {0, 1}) data[5 + i_dim] = con.stretch()[i_dim]; \
      Con_dir<Deformed_element> dir = con.direction(); \
      for (int i_side = 0; i_side < 2; ++i_side) { \
        data[7 + i_side] = dir.i_dim[i_side]; \
        data[9 + i_side] = dir.face_sign[i_side]; \
      } \
      if (con.order_reversed()) { \
        /* rearrange order to turn reversed connections into un-reversed ones, \
           since mesh file format does not support reversed connecions */ \
        for (int i : {7, 9}) std::swap(data[i], data[i + 1]); \
      } \
      h5_write_row(ref_con_dset, 11, start + i_con, data); \
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
}

void Accessible_mesh::read_file(std::string file_name) {
  H5::H5File file(file_name + ".mesh.h5", H5F_ACC_RDONLY);
  hsize_t dims [2];
  // read elements
  auto vert_pos_dset = file.openDataSet("/vertices/position");
  auto vert_ind_dset = file.openDataSet("/elements/vertices");
  auto nom_pos_dset = file.openDataSet("/elements/nominal_position");
  auto is_def_dset = file.openDataSet("/elements/is_deformed");
  auto ref_level_dset = file.openDataSet("/elements/refinement_level");
  auto aniso_ref_level_dset = file.openDataSet("/elements/aniso_refinement_level");
  is_def_dset.getSpace().getSimpleExtentDims(dims);
  int n_elem = dims[0];
  int n_vert = params.n_vertices();
  std::vector<Element*> elem_ptrs(n_elem); // really need to switch to storing a flat array of fully polymorphic elements to avoid this nonsense
  std::vector<Deformed_element*> def_elem_ptrs(n_elem, nullptr);
  for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
    std::vector<int> nom_pos(params.n_dim);
    h5_read_row(nom_pos_dset, params.n_dim, i_elem, nom_pos.data());
    std::vector<Int> np;
    for (int p : nom_pos) np.push_back(p);
    int ref_level = h5_read_value<int>(ref_level_dset, i_elem);
    int aniso_ref_level = h5_read_value<int>(aniso_ref_level_dset, i_elem);
    int is_def = h5_read_value<bool>(is_def_dset, i_elem);
    int sn = add_element(ref_level, is_def, np, tree ? tree->origin() : Mat<>::Zero(params.n_dim), aniso_ref_level);
    int vert_inds[8] {};
    h5_read_row(vert_ind_dset, n_vert, i_elem, vert_inds);
    auto& elem = element(ref_level, is_def, sn);
    elem_ptrs[i_elem] = &elem;
    if (is_def) def_elem_ptrs[i_elem] = &def.elems.at(ref_level, sn);
    for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
      Mat<3> p;
      h5_read_row(vert_pos_dset, params.n_dim, vert_inds[i_vert], p.data());
      elem.shape().vertex(i_vert).set_pos(p);
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
      _connect(def_el_ar, {{data[2], data[3]}, {bool(data[4]), bool(data[5])}});
      if (bool(def_el_ar[0]->tree) != bool(def_el_ar[1]->tree)) extrude_cons.push_back(def.cons.back().get());
    } else {
      _connect(el_ar, {data[2]});
    }
  }
  // read refined connections
  auto ref_con_dset = file.openDataSet("/connections/refined");
  ref_con_dset.getSpace().getSimpleExtentDims(dims);
  n_con = dims[0];
  for (int i_con = 0; i_con < n_con; ++i_con) {
    int data [11];
    h5_read_row(ref_con_dset, 11, i_con, data);
    std::array<bool, 2> stretch {bool(data[5]), bool(data[6])};
    int n_fine = math::pow(2, params.n_dim - 1 - stretch[0] - stretch[1]);
    Element* coarse = elem_ptrs[data[0]];
    std::vector<Element*> fine(n_fine);
    bool is_def = coarse->get_is_deformed();
    for (int i_fine = 0; i_fine < n_fine; ++i_fine) {
      fine[i_fine] = elem_ptrs[data[1 + i_fine]];
      is_def = is_def && fine[i_fine]->get_is_deformed();
    }
    Con_dir<Deformed_element> dir {{data[7], data[8]}, {bool(data[9]), bool(data[10])}};
    if (is_def) {
      Deformed_element* def_coarse = def_elem_ptrs[data[0]];
      std::vector<Deformed_element*> def_fine(n_fine);
      for (int i_fine = 0; i_fine < n_fine; ++i_fine) def_fine[i_fine] = def_elem_ptrs[data[1 + i_fine]];
      _connect(def_coarse, def_fine, dir, stretch);
    } else {
      _connect(coarse, fine, dir);
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
      def.bound_cons.emplace_back(new Typed_bound_connection<Deformed_element>(
        *def_elem_ptrs[data[0]], data[2], data[3], data[1], bound_conds[data[1]]->n_prescribed(params.n_dim)
      ));
    } else {
      car.bound_cons.emplace_back(new Typed_bound_connection<Element>(
        *elem_ptrs[data[0]], data[2], data[3], data[1], bound_conds[data[1]]->n_prescribed(params.n_dim)
      ));
    }
  }
  cleanup();
}

Accessible_mesh::Accessible_mesh(std::string file_name, std::vector<Flow_bc*> extremal_bcs, Turbulence_model turb,
                                 Surface_geom* geometry, Flow_bc* surface_bc)
: Accessible_mesh(read_params(file_name), read_root_sz(file_name), turb) {
  // take ownership of these to avoid memory leaks in case of exception
  std::unique_ptr<Flow_bc> fbc(surface_bc);
  std::unique_ptr<Surface_geom> g(geometry);
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
    surf_bc_sn = add_boundary_condition(fbc.release());
    surf_geom.reset(g.release());
  }
  read_file(file_name);
}

Accessible_mesh::Accessible_mesh(std::string file_name, std::vector<Flow_bc*> flow_bcs, Turbulence_model turb)
: Accessible_mesh(read_params(file_name), read_root_sz(file_name), turb) {
  for (unsigned i_bc = 0; i_bc < flow_bcs.size(); ++i_bc) add_boundary_condition(flow_bcs[i_bc]);
  read_file(file_name);
}

void write_polymesh_file(std::string dir_name, std::string name, std::string cls, int n_entries,
                         std::function<std::string(int)> entries, std::string note = "") {
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

void Accessible_mesh::export_polymesh(std::string dir_name) {
  dir_name = dir_name + "polyMesh/";
  if (std::filesystem::exists(dir_name)) std::filesystem::remove_all(dir_name);
  std::filesystem::create_directory(dir_name);
  auto verts = _blocks.verts();
  auto& bound_cons = boundary_connections();
  auto& elems = elements();
  auto& elem_cons = element_connections();
  int n_internal = elem_cons.size();
  int n_faces = n_internal + bound_cons.size();
  std::string face_note = format_str(200, "nPoints:%i nCells:%i nFaces:%i nInternalFaces:%i",
                                     verts.size(), elems.size(), n_faces, n_internal);
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) elems[i_elem].record = i_elem;
  write_polymesh_file(dir_name, "points", "vectorField", verts.size(), [&](int i_vert) {
    next::Vertex& vert = verts[i_vert];
    vert.record.clear();
    vert.record.push_back(i_vert);
    return format_str(100, "(%.20e %.20e %.20e)", vert.point({})[0], vert.point({})[1], vert.point({})[2]);
  });
  std::vector<int> owners(n_faces);
  std::vector<int> neighbors(n_internal);
  std::vector<std::vector<int>> faces(n_faces);
  int i_face = 0;
  auto add_verts = [&](Element& elem, int i_dim, int face_sign, bool flip) {
    auto& verts = faces[i_face++];
    for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
      if ((i_vert/math::pow(2, params.n_dim - i_dim - 1))%2 == face_sign) {
        verts.push_back(elem.shape().vertex(i_vert).record[0]);
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

void Accessible_mesh::visualize(std::string format, std::string file_name, double time) {
  next::Sequence<const next::Block&> shapes(
    [this](std::size_t index)->const next::Block& {return elements()[index].shape();},
    [this]()->std::size_t {return elements().size();}
  );
  next::Block::visualize(format, file_name, shapes, time);
}

}
