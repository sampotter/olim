#include "moore_marcher.hpp"

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
