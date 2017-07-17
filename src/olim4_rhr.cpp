#include "olim4_rhr.hpp"

#include "olim_util.hpp"

void olim4_rhr::update_node_value_impl(int i, int j, double & T) {
  node * nb[4] = {0x0, 0x0, 0x0, 0x0}; // NESW
  get_valid_neighbors(i, j, (abstract_node **) nb);
  double s = S(i, j), h = get_h(), sh = s*h;
  node * x0 = 0x0;
  node * x1 = 0x0;

  for (int k = 0; k < 4; ++k) {
    if ((x0 = nb[k]) && !nb[(k + 1) % 4] && !nb[(k + 3) % 4]) {
      T = std::min(T, x0->get_value() + sh);
    }
  }

  for (int k = 0; k < 4; ++k) {
    if ((x0 = nb[k]) && (x1 = nb[(k + 1) % 4])) {
      T = std::min(T, rhr_adj(x0->get_value(), x1->get_value(), s, h));
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
