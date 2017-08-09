#ifndef __OLIM8_IMPL_HPP__
#define __OLIM8_IMPL_HPP__

#include <algorithm>

template <class update_rules>
void olim8<update_rules>::update_impl(int i, int j, double & T) {
  abstract_node * nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  abstract_node * x0 = 0x0;
  abstract_node * x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double s = S(i, j), h = get_h();

  /*
   * First, do the adjacent and diagonal single point updates (and
   * only if the update triangles upon which the update edges are
   * incident are *not* present):
   */
  for (int k = 0; k < 8; k += 2) { // adjacent
    if ((x0 = nb[k]) && !nb[(k + 6) % 8] && !nb[(k + 7) % 8] &&
        !nb[(k + 1) % 8] && !nb[(k + 2) % 8]) {
      T = std::min(T, this->adj1pt(x0->get_value(), s, h));
    }
  }
  for (int k = 1; k < 8; k += 2) { // diagonal
    if ((x0 = nb[k]) && !nb[(k + 7) % 8] && !nb[(k + 1) % 8]) {
      T = std::min(T, this->diag1pt(x0->get_value(), s, h));
    }
  }

  /**
   * Next, do the diagonal triangle updates.
   */
  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && (x1 = nb[k + 1])) {
      T = std::min(T, this->diag2pt(x0->get_value(), x1->get_value(), s, h));
    }
  }
  for (int k = 1; k < 8; k += 2) {
    if ((x0 = nb[(k + 1) % 8]) && (x1 = nb[k])) {
      T = std::min(T, this->diag2pt(x0->get_value(), x1->get_value(), s, h));
    }
  }

  /**
   * Finally, do the adjacent triangle updates.
   */
  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && (x1 = nb[(k + 2) % 8])) {
      T = std::min(T, this->adj2pt(x0->get_value(), x1->get_value(), s, h));
    }
  }
}

#endif // __OLIM8_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
