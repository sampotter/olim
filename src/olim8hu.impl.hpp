#ifndef __OLIM8HU_IMPL_HPP__
#define __OLIM8HU_IMPL_HPP__

#include <algorithm>
#include <cassert>

template <class update_rules>
void olim8hu<update_rules>::update_impl(int i, int j, double & T) {
  abstract_node * nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  abstract_node * x0 = 0x0;
  abstract_node * x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double s = speed(i, j), h = get_h(), u0, Tnew;

  /**
   * Do every one-point update, and keep track of the index of the
   * minimizing point.
   */
  int argmin = -1;
  for (int k = 0; k < 8; ++k) {
    if ((x0 = nb[k])) {
      u0 = x0->get_value();
      Tnew = k % 2 == 0 ? this->adj1pt(u0, s, h) : this->diag1pt(u0, s, h);
      if (Tnew < T) {
        T = Tnew;
        argmin = k;
      }
    }
  }
  assert(argmin != -1);

  /**
   * Do the two-point updates for the triangles which are adjacent to
   * the point that had the minimum one-point update.
   */
  int kprev = (argmin + 7) % 8, knext = (argmin + 1) % 8;
  u0 = nb[argmin]->get_value();
  if (argmin % 2 == 0) {
    if ((x1 = nb[kprev])) {
      T = std::min(T, this->diag2pt(u0, x1->get_value(), s, h));
    }
    if ((x1 = nb[knext])) {
      T = std::min(T, this->diag2pt(u0, x1->get_value(), s, h));
    }
  } else {
    x1 = nb[kprev];
    abstract_node * x2 = nb[knext];
    if (x1) {
      T = std::min(T, this->diag2pt(x1->get_value(), u0, s, h));
    }
    if (x2) {
      T = std::min(T, this->diag2pt(x2->get_value(), u0, s, h));
    }
    if (x1 && x2) {
      T = std::min(T, this->adj2pt(x1->get_value(), x2->get_value(), s, h));
    }
  }
}

#endif // __OLIM8HU_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
