#include "olim8_mp0.hpp"

#include <cmath>

#include "olim_util.hpp"

void olim8_mp0::update_impl(int i, int j, double & T) {
  abstract_node * nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  abstract_node * x0 = 0x0;
  abstract_node * x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double h = get_h(), u0, u1, s = S(i, j), s_est;

  /*
   * First, do the adjacent and diagonal single point updates (and
   * only if the update triangles upon which the update edges are
   * incident are *not* present):
   */
  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 6) % 8] && !nb[(k + 7) % 8] &&
        !nb[(k + 1) % 8] && !nb[(k + 2) % 8]) {
      s_est = (s + S(i + di[k], j + dj[k]))/2;
      T = std::min(T, x0->get_value() + h*s_est);
    }
  }
  for (int k = 1; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 7) % 8] && !nb[(k + 1) % 8]) {
      s_est = (s + S(i + di[k], j + dj[k]))/2;
      T = std::min(T, x0->get_value() + h*s_est*std::sqrt(2));
    }
  }

  /**
   * Next, do the diagonal triangle updates.
   */
  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && (x1 = nb[k + 1])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s_est = get_s_est(s, i + di[k], j + dj[k], i + di[k + 1], j + dj[k + 1]);
      T = std::min(T, rhr_diag(u0, u1, s_est, h));
    }
  }
  for (int k = 1, l = k + 1; k < 8; k += 2, l = (k + 1) % 8) {
    if ((x0 = nb[l]) && (x1 = nb[k])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s_est = get_s_est(s, i + di[l], j + dj[l], i + di[k], j + dj[k]);
      T = std::min(T, rhr_diag(u0, u1, s_est, h));
    }
  }

  /**
   * Finally, do the adjacent triangle updates.
   */
  for (int k = 0, l = k + 2; k < 8; k += 2, l = (k + 2) % 8) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s_est = get_s_est(s, i + di[k], j + dj[k], i + di[l], j + dj[l]);
      T = std::min(T, rhr_adj(u0, u1, s_est, h));
    }
  }
}

double olim8_mp0::get_s_est(double s, int i0, int j0, int i1, int j1) {
  return (s + (S(i0, j0) + S(i1, j1))/2)/2;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
