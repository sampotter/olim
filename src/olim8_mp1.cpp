#include "olim8_mp1.hpp"

#include <cmath>

#include "olim8_util.hpp"

void olim8_mp1::update_node_value_impl(int i, int j, double & T) {
  node* nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  node *x0 = 0x0, *x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double h = get_h(), s = S(i, j), u0, u1, sbar0, sbar1;

  /*
   * First, do the adjacent and diagonal single point updates (and
   * only if the update triangles upon which the update edges are
   * incident are *not* present):
   */
  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 6) % 8] && !nb[(k + 7) % 8] &&
        !nb[(k + 1) % 8] && !nb[(k + 2) % 8]) {
      sbar0 = (s + S(i + di[k], j + dj[k]))/2;
      T = std::min(T, x0->get_value() + h*sbar0);
    }
  }
  for (int k = 1; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 7) % 8] && !nb[(k + 1) % 8]) {
      sbar0 = (s + S(i + di[k], j + dj[k]))/2;
      T = std::min(T, x0->get_value() + h*sbar0*std::sqrt(2));
    }
  }

  /**
   * Next, do the diagonal triangle updates.
   */
  for (int k = 0, l = k + 1; k < 8; k += 2, l = k + 1) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      sbar0 = (s + S(i + di[k], j + dj[k]))/2;
      sbar1 = (s + S(i + di[l], j + dj[l]))/2;
      T = std::min(T, sbar0 == sbar1 ? rhr_diag(u0, u1, sbar0, h) :
                   mp1_diag(u0, u1, sbar0, sbar1, h));
    }
  }
  for (int k = 1, l = k + 1; k < 8; k += 2, l = (k + 1) % 8) {
    if ((x0 = nb[l]) && (x1 = nb[k])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      sbar0 = (s + S(i + di[l], j + dj[l]))/2;
      sbar1 = (s + S(i + di[k], j + dj[k]))/2;
      T = std::min(T, sbar0 == sbar1 ? rhr_diag(u0, u1, sbar0, h) :
                   mp1_diag(u0, u1, sbar0, sbar1, h));
    }
  }

  /**
   * Finally, do the adjacent triangle updates.
   */
  for (int k = 0, l = k + 2; k < 8; k += 2, l = (k + 2) % 8) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      sbar0 = (s + S(i + di[k], j + dj[k]))/2;
      sbar1 = (s + S(i + di[l], j + dj[l]))/2;
      T = std::min(T, sbar0 == sbar1 ? rhr_adj(u0, u1, sbar0, h) :
                   mp1_adj(u0, u1, sbar0, sbar1, h));
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
