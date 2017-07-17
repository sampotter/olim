#include "olim4_mp0.hpp"

#include "olim_util.hpp"

void olim4_mp0::update_node_value_impl(int i, int j, double & T) {
  node * nb[4] = {0x0, 0x0, 0x0, 0x0}; // NESW
  get_valid_neighbors(i, j, (abstract_node **) nb);
  double s = S(i, j), h = get_h(), s_est;
  node * x0 = 0x0;
  node * x1 = 0x0;

  for (int k = 0; k < 4; ++k) {
    if ((x0 = nb[k]) && !nb[(k + 1) % 4] && !nb[(k + 3) % 4]) {
      s_est = (s + S(i + di[k], j + dj[k]))/2;
      T = std::min(T, x0->get_value() + h*s_est);
    }
  }

  for (int k = 0, l = 1; k < 4; ++k, l = (k + 1) % 4) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      s_est = (s + (S(i + di[k], j + dj[k]) + S(i + di[l], j + dj[l]))/2)/2;
      T = std::min(T, rhr_adj(x0->get_value(), x1->get_value(), s_est, h));
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
