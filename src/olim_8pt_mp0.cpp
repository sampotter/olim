#include "olim_8pt_mp0.hpp"

#include <cmath>

#include "olim_8pt_util.hpp"

void olim_8pt_mp0::update_node_value_impl(size_t i, size_t j, double & T) {
  static int di[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
  static int dj[8] = {0, 1, 1, 1, 0, -1, -1, -1};
  node* nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  node *x0 = 0x0, *x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double h = get_h(), x = h*i, y = h*j, u0, u1;
  double s = S(x, y), s_est;

  /*
   * First, do the adjacent and diagonal single point updates (and
   * only if the update triangles upon which the update edges are
   * incident are *not* present):
   */
  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 6) % 8] && !nb[(k + 7) % 8] &&
        !nb[(k + 1) % 8] && !nb[(k + 2) % 8]) {
      s_est = (s + S(x + h*di[k], x + h*dj[k]))/2;
      T = std::min(T, x0->get_value() + h*s_est);
    }
  }
  for (int k = 1; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 7) % 8] && !nb[(k + 1) % 8]) {
      s_est = (s + S(x + h*di[k], x + h*dj[k]))/2;
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
      s_est = get_s_est(s, x + h*di[k], y + h*dj[k], h*di[k + 1], h*dj[k + 1]);
      T = std::min(T, rhr_diag(u0, u1, s_est, h));
    }
  }
  for (int k = 1; k < 8; k += 2) {
    if ((x0 = nb[(k + 1) % 8]) && (x1 = nb[k])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s_est = get_s_est(s, x + h*di[k + 1], y + h*dj[k + 1], h*di[k], h*dj[k]);
      T = std::min(T, rhr_diag(u0, u1, s_est, h));
    }
  }

  /**
   * Finally, do the adjacent triangle updates.
   */
  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && (x1 = nb[(k + 2) % 8])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s_est = get_s_est(s, x + h*di[k], y + h*dj[k], h*di[k + 2], h*dj[k + 2]);
      T = std::min(T, rhr_adj(u0, u1, s_est, h));
    }
  }
}

double olim_8pt_mp0::get_s_est(double s, double x0, double y0, double x1,
                               double y1) const {
  return (s + (S(x0, y0) + S(x1, y1))/2)/2;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
