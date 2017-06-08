#include "moore_marcher.hpp"

void moore_marcher::stage_neighbors_impl(size_t i, size_t j) {
  static int offsets[8][2] = {{-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1},
                              {0, -1}, {-1, -1}};

  for (int k = 0; k < 8; ++k) {
    stage_neighbor(i + offsets[k][0], j + offsets[k][1]);
  }

  size_t a, b;
  for (int k = 0; k < 8; ++k) {
    a = i + offsets[k][0], b = j + offsets[k][1];
    if (in_bounds(a, b) && !this->operator()(a, b).is_valid()) {
      update_node_value(a, b);
    }
  }
}

void moore_marcher::get_valid_neighbors(size_t i, size_t j, node ** nb) {
  if (is_valid(i - 1, j)) nb[0] = &this->operator()(i - 1, j);         // N
  if (is_valid(i - 1, j + 1)) nb[0] = &this->operator()(i - 1, j + 1); // NE
  if (is_valid(i, j + 1)) nb[0] = &this->operator()(i, j + 1);         // E
  if (is_valid(i + 1, j + 1)) nb[0] = &this->operator()(i + 1, j + 1); // SE
  if (is_valid(i + 1, j)) nb[0] = &this->operator()(i + 1, j);         // S
  if (is_valid(i + 1, j - 1)) nb[0] = &this->operator()(i + 1, j - 1); // SW
  if (is_valid(i, j - 1)) nb[0] = &this->operator()(i, j - 1);         // W
  if (is_valid(i - 1, j - 1)) nb[0] = &this->operator()(i - 1, j - 1); // NW
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
