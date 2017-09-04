#include "olim4_mp0.hpp"

#include "olim_util.hpp"

void olim4_mp0::update_impl(int i, int j, double & T) {
  abstract_node * nb[4] = {0x0, 0x0, 0x0, 0x0}; // NESW
  get_valid_neighbors(i, j, nb);
  double s = speed(i, j), h = get_h(), s_est;

  for (int k = 0; k < 4; ++k) {
    if (nb[k] && !nb[(k + 1) % 4] && !nb[(k + 3) % 4]) {
      s_est = (s + speed(i + di[k], j + dj[k]))/2;
      T = std::min(T, nb[k]->get_value() + h*s_est);
    }
  }

  for (int k = 0, l = 1; k < 4; ++k, l = (k + 1) % 4) {
    if (nb[k] && nb[l]) {
      s_est = (s + (speed(i + di[k], j + dj[k]) + speed(i + di[l], j + dj[l]))/2)/2;
      T = std::min(T, rhr_adj(nb[k]->get_value(), nb[l]->get_value(), s_est, h));
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
