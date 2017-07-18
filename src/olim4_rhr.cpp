#include "olim4_rhr.hpp"

#include "olim_util.hpp"

void olim4_rhr::update_node_value_impl(int i, int j, double & T) {
  abstract_node * nb[4] = {0x0, 0x0, 0x0, 0x0}; // NESW
  get_valid_neighbors(i, j, nb);
  double s = S(i, j), h = get_h(), sh = s*h;

  for (int k = 0; k < 4; ++k) {
    if (nb[k] && !nb[(k + 1) % 4] && !nb[(k + 3) % 4]) {
      T = std::min(T, nb[k]->get_value() + sh);
    }
  }

  for (int k = 0, l = k + 1; k < 4; ++k, l = (k + 1) % 4) {
    if (nb[k] && nb[l]) {
      T = std::min(T, rhr_adj(nb[k]->get_value(), nb[l]->get_value(), s, h));
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
