#include "stats.hpp"

#include <algorithm>
#include <cassert>

olim3d_node_stats::olim3d_node_stats() {
  std::fill(
    _num_line_updates,
    _num_line_updates + static_cast<int>(line_update::NUM),
    0);
  std::fill(
    _num_tri_updates,
    _num_tri_updates + static_cast<int>(tri_update::NUM),
    0);
  std::fill(
    _num_degenerate_tri_updates,
    _num_degenerate_tri_updates + static_cast<int>(tri_update::NUM),
    0);
  std::fill(
    _num_tetra_updates,
    _num_tetra_updates + static_cast<int>(tetra_update::NUM),
    0);
}

int olim3d_node_stats::num_line_updates() const {
  return _num_line_updates[0] + _num_line_updates[1] + _num_line_updates[2];
}

int olim3d_node_stats::num_line_updates(int d) const {
  assert(d == 1 || d == 2 || d == 3);
  return _num_line_updates[d - 1];
}

void olim3d_node_stats::inc_line_updates(int d) {
  assert(d == 1 || d == 2 || d == 3);
  ++_num_line_updates[d - 1];
}

int olim3d_node_stats::num_tri_updates() const {
  int tmp = 0;
  for (int i = 0; i < static_cast<int>(tri_update::NUM); ++i) {
    tmp += _num_tri_updates[i];
  }
  return tmp;
}

int olim3d_node_stats::num_degenerate_tri_updates() const {
  int tmp = 0;
  for (int i = 0; i < static_cast<int>(tri_update::NUM); ++i) {
    tmp += _num_degenerate_tri_updates[i];
  }
  return tmp;
}

static int get_tri_index(int d1, int d2) {
  int i = -1;
  if (d1 > d2) std::swap(d1, d2);
  if (d1 == 1) {
    if (d2 == 1) i = 0;
    else if (d2 == 2) i = 1;
    else if (d2 == 3) i = 2;
    else assert(false);
  } else if (d1 == 2) {
    if (d2 == 2) i = 3;
    else if (d2 == 3) i = 4;
    else assert(false);
  } else {
    assert(false);
  }
  return i;
}

void olim3d_node_stats::inc_tri_updates(int d1, int d2, bool degenerate) {
  int i = get_tri_index(d1, d2);
  if (degenerate) ++_num_degenerate_tri_updates[i];
  ++_num_tri_updates[i];
}

int olim3d_node_stats::num_tetra_updates() const {
  int tmp = 0;
  for (int i = 0; i < static_cast<int>(tetra_update::NUM); ++i) {
    tmp += _num_tetra_updates[i];
  }
  return tmp;
}

void olim3d_node_stats::inc_tetra_updates(int d1, int d2, int d3) {
  int a = std::min(d1, std::min(d2, d3));
  int b = std::max(d1, std::max(d2, d3));
  d2 += d1 + d3 - a - b;
  d1 = a;
  d3 = b;
  if (d1 == 1) {
    if (d2 == 1) {
      if (d3 == 1) ++_num_tetra_updates[0];
      else if (d3 == 2) ++_num_tetra_updates[1];
      else if (d3 == 3) ++_num_tetra_updates[2];
      else assert(false);
    } else if (d2 == 2) {
      if (d3 == 2) ++_num_tetra_updates[3];
      else if (d3 == 3) ++_num_tetra_updates[4];
      else assert(false);
    } else {
      assert(false);
    }
  } else if (d1 == 2) {
    if (d2 == 2) {
      if (d3 == 2) ++_num_tetra_updates[5];
      else if (d3 == 3) ++_num_tetra_updates[6];
      else assert(false);
    } else {
      assert(false);
    }
  } else {
    assert(false);
  }
}
