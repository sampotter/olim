#include "moore_marcher.hpp"

// neighbor order: N, NE, E, SE, S, SW, W, NW (clockwise from north)

int moore_marcher::di[] = {-1, -1, 0, 1, 1, 1, 0, -1};
int moore_marcher::dj[] = {0, 1, 1, 1, 0, -1, -1, -1};

void moore_marcher::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();

  for (int k = 0; k < 8; ++k) {
    stage_neighbor(i + di[k], j + dj[k]);
  }

  int a, b;
  for (int k = 0; k < 8; ++k) {
    a = i + di[k];
    b = j + dj[k];
    if (in_bounds(a, b) && !this->operator()(a, b).is_valid()) {
      update_node_value(a, b);
    }
  }
}

void moore_marcher::get_valid_neighbors(int i, int j, node ** nb) {
  int a, b;
  for (int k = 0; k < 8; ++k) {
    a = i + di[k];
    b = j + dj[k];
    if (is_valid(a, b)) {
      nb[k] = &this->operator()(a, b);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
