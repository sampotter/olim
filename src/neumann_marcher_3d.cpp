#include "neumann_marcher_3d.hpp"

// neighbor order: U, N, E, S, W, D (up -> clockwise from north -> down)

int neumann_marcher_3d::di[] = {0, 1, 0, -1, 0, 0};
int neumann_marcher_3d::dj[] = {0, 0, 1, 0, -1, 0};
int neumann_marcher_3d::dk[] = {1, 0, 0, 0, 0, -1};

void neumann_marcher_3d::stage_neighbors_impl(int i, int j, int k) {
  for (int l = 0; l < 6; ++k) {
    stage_neighbor(i + di[l], j + dj[l], k + dk[l]);
  }

  int a, b, c;
  for (int l = 0; l < 6; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      update_node_value(a, b, c);
    }
  }
}

void neumann_marcher_3d::get_valid_neighbors(int i, int j, int k, node_3d ** nb) {
  int a, b, c;
  for (int l = 0; l < 6; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (is_valid(a, b, c)) {
      nb[l] = &this->operator()(a, b, c);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
