#include "olim8_mp0l.hpp"

#include <cmath>

#include "olim_util.hpp"

void olim8_mp0l::update_impl(int i, int j, double & T) {
  abstract_node * nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  abstract_node * x0 = 0x0;
  abstract_node * x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double h = get_h(), u0, u1, s = S(i, j), s0, s1;

  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 6) % 8] && !nb[(k + 7) % 8] &&
        !nb[(k + 1) % 8] && !nb[(k + 2) % 8]) {
      s0 = S(i + di[k], j + dj[k]);
      T = std::min(T, x0->get_value() + h*(s + s0)/2);
    }
  }
  for (int k = 1; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 7) % 8] && !nb[(k + 1) % 8]) {
      s0 = S(i + di[k], j + dj[k]);
      T = std::min(T, x0->get_value() + h*(s + s0)*std::sqrt(2)/2);
    }
  }

  for (int k = 0, l = 1; k < 8; k += 2, l = k + 1) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s0 = S(i + di[k], j + dj[k]);
      s1 = S(i + di[l], j + dj[l]);
      T = std::min(T, mp0l_diag(u0, u1, s, s0, s1, h));
    }
  }
  for (int k = 1, l = k + 1; k < 8; k += 2, l = (k + 1) % 8) {
    if ((x0 = nb[l]) && (x1 = nb[k])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s0 = S(i + di[l], j + dj[l]);
      s1 = S(i + di[k], j + dj[k]);
      T = std::min(T, mp0l_diag(u0, u1, s, s0, s1, h));
    }
  }

  for (int k = 0, l = k + 2; k < 8; k += 2, l = (k + 2) % 8) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s0 = S(i + di[k], j + dj[k]);
      s1 = S(i + di[l], j + dj[l]);
      T = std::min(T, mp0l_adj(u0, u1, s, s0, s1, h));
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
