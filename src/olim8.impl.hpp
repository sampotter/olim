#ifndef __OLIM8_IMPL_HPP__
#define __OLIM8_IMPL_HPP__

#include <algorithm>

#include <src/config.hpp>

template <class update_rules>
void olim8<update_rules>::update_impl(int i, int j, double & T) {
  abstract_node * nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  abstract_node * x0 = 0x0;
  abstract_node * x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double s = speed(i, j), s0, s1, u0, u1, h = get_h();

  /*
   * First, do the adjacent and diagonal single point updates (and
   * only if the update triangles upon which the update edges are
   * incident are *not* present):
   */
  for (int k = 0; k < 8; k += 2) { // adjacent
    if ((x0 = nb[k]) && !nb[(k + 6) % 8] && !nb[(k + 7) % 8] &&
        !nb[(k + 1) % 8] && !nb[(k + 2) % 8]) {
      u0 = x0->get_value();
      s0 = speed(i + di[k], j + dj[k]);
      T = std::min(T, this->template line<1>(u0, s, s0, h));
    }
  }
  for (int k = 1; k < 8; k += 2) { // diagonal
    if ((x0 = nb[k]) && !nb[(k + 7) % 8] && !nb[(k + 1) % 8]) {
      u0 = x0->get_value();
      s0 = speed(i + di[k], j + dj[k]);
      T = std::min(T, this->template line<2>(u0, s, s0, h));
    }
  }

  /**
   * Next, do the diagonal triangle updates.
   */
  for (int k = 0, l = 1; k < 8; k += 2, l = k + 1) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s0 = speed(i + di[k], j + dj[k]);
      s1 = speed(i + di[l], j + dj[l]);
      T = std::min(T, this->tri12(u0, u1, s, s0, s1, h));
    }
  }
  for (int k = 1, l = 2; k < 8; k += 2, l = (k + 1) % 8) {
    if ((x0 = nb[l]) && (x1 = nb[k])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s0 = speed(i + di[l], j + dj[l]);
      s1 = speed(i + di[k], j + dj[k]);
      T = std::min(T, this->tri12(u0, u1, s, s0, s1, h));
    }
  }

#ifdef OLIM8_ADJ_UPDATES
  /**
   * Finally, do the adjacent triangle updates.
   */
  for (int k = 0, l = 2; k < 8; k += 2, l = (k + 2) % 8) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s0 = speed(i + di[k], j + dj[k]);
      s1 = speed(i + di[l], j + dj[l]);
      T = std::min(T, this->tri11(u0, u1, s, s0, s1, h));
    }
  }
#endif // OLIM8_ADJ_UPDATES
}

#endif // __OLIM8_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
