#include "basic_marcher.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

void basic_marcher::update_node_value_impl(int i, int j, double & T) {
  abstract_node * nb[4] = {nullptr, nullptr, nullptr, nullptr};
  get_valid_neighbors(i, j, nb);
  double sh = get_h()*S(i, j);
  double T1 = 0, T2 = 0, disc = 0;
  for (int k = 0, k1 = 1; k < 4; ++k, k1 = (k1 + 1) % 4) {
    if (nb[k] && nb[k1]) {
      T1 = nb[k]->get_value();
      T2 = nb[k1]->get_value();
      disc = 2*sh*sh - (T1 - T2)*(T1 - T2);
      T = disc > 0 ? std::min(T, (T1 + T2 + std::sqrt(disc))/2) : T;
    } else if (nb[k]) {
      T = std::min(T, nb[k]->get_value() + sh);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
