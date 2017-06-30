#include "neumann_marcher.hpp"

void neumann_marcher::stage_neighbors_impl(int i, int j) {
  static int offsets[4][2] = {{-1, 0}, {0, -1}, {0, 1}, {1, 0}};

  for (int k = 0; k < 4; ++k) {
    stage_neighbor(i + offsets[k][0], j + offsets[k][1]);
  }

  int a, b;
  for (int k = 0; k < 4; ++k) {
    a = i + offsets[k][0], b = j + offsets[k][1];
    if (in_bounds(a, b) && !this->operator()(a, b).is_valid()) {
      update_node_value(a, b);
    }
  }
}

void neumann_marcher::get_valid_neighbors(int i, int j, node ** nb) {
  if (is_valid(i - 1, j)) nb[0] = &this->operator()(i - 1, j); // N
  if (is_valid(i, j + 1)) nb[1] = &this->operator()(i, j + 1); // E
  if (is_valid(i + 1, j)) nb[2] = &this->operator()(i + 1, j); // S
  if (is_valid(i, j - 1)) nb[3] = &this->operator()(i, j - 1); // W
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
