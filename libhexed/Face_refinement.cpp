#include <hexed/Face_refinement.hpp>
#include <hexed/vertex_inds.hpp>

namespace hexed {

Face_refinement::Face_refinement(Face& face, int sdim)
: _coarse{this}
, _split_dim{sdim}
, _fine0{face.storage_params(), face.i_dim(), face.sign(), face.is_deformed()}
, _fine1{face.storage_params(), face.i_dim(), face.sign(), face.is_deformed()}
{
  HEXED_ASSERT(0 <= sdim && sdim < face.storage_params().n_dim - 1, "`split_dim` out of bounds")
  face.connect(_coarse);
  _fine0.associate(*this);
  _fine1.associate(*this);
}

Face& find_element_face(Face& face, bool upstream) {
  if (upstream) {
    HEXED_ASSERT(face.associated(), "Face is not associated.")
    if (face.element()) {
      return face;
    } else {
      return find_element_face(face.face_ref_coarse()->coarse(), true);
    }
  } else {
    HEXED_ASSERT(face.connected(), "Face is not connected.")
    if (face.neighbor_connection()) {
      return find_element_face(face.neighbor_connection()->opposite_face(face), true);
    } else {
      return find_element_face(*face.face_ref_fine()->fine()[0], false);
    }
  }
}

void find_elements(Face& face, std::vector<Element*>& elems, std::vector<int> bounds, bool upstream) {
  if (upstream) {
    HEXED_ASSERT(face.associated(), "Face is not associated.")
    if (face.element()) {
      for (int i = bounds[0]; i < bounds[1]; ++i) {
        for (int j = bounds[2]; j < bounds[3]; ++j) {
          elems[2*i + j] = face.element();
        }
      }
    } else {
      find_elements(face.face_ref_coarse()->coarse(), elems, bounds, true);
    }
  } else {
    HEXED_ASSERT(face.connected(), "Face is not connected.")
    if (face.neighbor_connection()) {
      find_elements(face.neighbor_connection()->opposite_face(face), elems, bounds, true);
    } else {
      for (int i_fine = 0; i_fine < 2; ++i_fine) {
        auto& ref = *face.face_ref_fine();
        auto b = bounds;
        b[2*(ref.split_dim() + (ref.coarse().storage_params().n_dim < 3)) + !i_fine] += math::sign(i_fine);
        find_elements(*ref.fine()[i_fine], elems, b, false);
      }
    }
  }
}

std::array<std::vector<Element*>, 2> Face_refinement::elements() {
  std::array<std::vector<Element*>, 2> elems;
  std::array<Face*, 2> start_from {fine()[0], &coarse()};
  auto par = coarse().storage_params();
  int n_fine = par.n_vertices()/2;
  auto dir_rev = _dir_reverse();
  for (int upstream = 0; upstream < 2; ++upstream) {
    int i_side = upstream != dir_rev.second;
    elems[i_side].resize(n_fine, nullptr);
    auto& elem_face = find_element_face(*start_from[upstream], upstream);
    std::vector<int> bounds {0, 1 + (par.n_dim > 2), 0, 1 + (par.n_dim > 1)};
    find_elements(elem_face, elems[i_side], bounds, false);
  }
  std::vector<int> fvi = face_vertex_inds(par.n_dim, dir_rev.first);
  std::array<std::vector<Element*>, 2> permuted;
  for (int i_side = 0; i_side < 2; ++i_side) {
    permuted[i_side].resize(n_fine, nullptr);
    for (int i_fine = 0; i_fine < n_fine; ++i_fine) {
      if (dir_rev.second) permuted[i_side][i_fine] = elems[i_side][fvi[i_fine]];
      else permuted[i_side][fvi[i_fine]] = elems[i_side][i_fine];
    }
    for (Element* p : permuted[i_side]) HEXED_ASSERT(p, "null element")
  }
  return permuted;
}

Connection_direction Face_refinement::get_direction() {
  return _dir_reverse().first;
}

std::pair<Connection_direction, bool> Face_refinement::_dir_reverse() {
  HEXED_ASSERT(_fine0.connected(), "Not connected.")
  if (_fine0.neighbor_connection()) {
    return {_fine0.neighbor_connection()->get_direction(), &_fine0.neighbor_connection()->face(1) == &_fine0};
  }
  return _fine0.face_ref_fine()->_dir_reverse();
}

Kernel_face_refinement Face_refinement::kernel_face_refinement() {
  Kernel_face_refinement kfr;
  for (int is_ldg = 0; is_ldg < 2; ++is_ldg) {
    kfr.coarse[is_ldg] = coarse().flow_state()(is_ldg).data();
    for (int i_fine = 0; i_fine < 2; ++i_fine) {
      kfr.fine[i_fine][is_ldg] = fine()[i_fine]->flow_state()(is_ldg).data();
    }
  }
  kfr.mask = coarse().mask();
  kfr.split_dim = split_dim();
  return kfr;
}

}
