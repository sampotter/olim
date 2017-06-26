#include "basic_marcher.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

void basic_marcher::update_node_value_impl(size_t i, size_t j, double & T) {
  node* n = 0x0;
  node* nb[4] = {0x0, 0x0, 0x0, 0x0}; // NESW
  get_valid_neighbors(i, j, nb);
  double sh = get_h()*S(i, j);
  double tmp = 0, T1 = 0, T2 = 0, disc = 0;
  for (int k = 0, k1 = 1; k < 4; ++k, k1 = (k1 + 1) % 4) {
    if (nb[k] && nb[k1]) {
      T1 = nb[k]->get_value();
      T2 = nb[k1]->get_value();
      disc = 2*sh*sh - (T1 - T2)*(T1 - T2);
      if (disc <= 0) continue;
      tmp = (T1 + T2 + std::sqrt(disc))/2;
      if (tmp >= T1 && tmp >= T2) T = std::min(T, tmp);
    } else if (nb[k] || nb[k1]) {
      n = nb[k] ? nb[k] : nb[k1];
      T = std::min(T, n->get_value() + sh);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
