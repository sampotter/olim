#include "basic_marcher.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

void basic_marcher::update_node_value_impl(size_t i, size_t j) {
  node* n = nullptr;
  node* nb[4] = {nullptr, nullptr, nullptr, nullptr}; // NESW
  get_valid_neighbors(i, j, nb);
  double rhs = get_h()/F(i, j);
  double T = std::numeric_limits<double>::infinity();
  double tmp = 0, T1 = 0, T2 = 0, disc = 0;

  for (int k = 0, k1 = 1; k < 4; ++k, k1 = (k1 + 1) % 4) {
    if (nb[k] != nullptr && nb[k1] != nullptr) {
      T1 = nb[k]->get_value();
      T2 = nb[k1]->get_value();
      disc = 2*std::pow(rhs, 2) - std::pow(T1 - T2, 2);
      if (disc <= 0) continue;
      tmp = (T1 + T2 + std::sqrt(disc))/2;
      if (tmp >= T1 && tmp >= T2) T = std::min(T, tmp);
    } else if (nb[k] != nullptr || nb[k1] != nullptr) {
      n = nb[k] != nullptr ? nb[k] : nb[k1];
      T = std::min(T, n->get_value() + rhs);
    }
  }

  n = &this->operator()(i, j);
  assert(n->is_trial());
  n->set_value(T);
  adjust_heap_entry(n);
}

void basic_marcher::update_neighbors_impl(size_t i, size_t j) {
  static int offsets[4][2] = {{-1, 0}, {0, -1}, {0, 1}, {1, 0}};
  size_t a, b;

  // Change all 'far' neighbors to 'trial'
  node * n = nullptr;
  for (int k = 0; k < 4; ++k) {
    a = i + offsets[k][0], b = j + offsets[k][1];
    if (valid_index(a, b) && this->operator()(a, b).is_far()) {
      n = &this->operator()(a, b);
      n->set_trial();
      insert_into_heap(n);
    }
  }

  for (int k = 0; k < 4; ++k) {
    a = i + offsets[k][0], b = j + offsets[k][1];
    if (valid_index(a, b) && !this->operator()(a, b).is_valid()) {
      update_node_value(a, b);
    }
  }
}

void basic_marcher::get_valid_neighbors(size_t i, size_t j, node ** nb) {
  auto const is_good = [this] (size_t a, size_t b) {
    return valid_index(a, b) && this->operator()(a, b).is_valid();
  };
  if (is_good(i - 1, j)) nb[0] = &this->operator()(i - 1, j); // north ...
  if (is_good(i, j + 1)) nb[1] = &this->operator()(i, j + 1); // east ...
  if (is_good(i + 1, j)) nb[2] = &this->operator()(i + 1, j); // south ...
  if (is_good(i, j - 1)) nb[3] = &this->operator()(i, j - 1); // west ...
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
