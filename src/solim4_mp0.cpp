#include "solim4_mp0.hpp"

#include "olim_util.hpp"

void solim4_mp0::update_impl(int i, int j, double & T) {
  smart_node * n = &operator()(i, j);
  abstract_node * nb[4] = {0x0, 0x0, 0x0, 0x0}; // NESW
  get_valid_neighbors(i, j, nb);
  double s = S(i, j), h = get_h(), s_est, T_new, lam;

  for (int k = 0; k < 4; ++k) {
    if (!nb[k] || nb[(k + 1) % 4] || nb[(k + 3) % 4]) {
      continue;
    }
    s_est = (s + S(i + di[k], j + dj[k]))/2;
    T_new = nb[k]->get_value() + h*s_est;
    if (T_new < T) {
      T = T_new;
      n->set_parent_nodes({nb[k], nullptr});
      n->set_lambda(-1);
    }
  }

  for (int k = 0, l = 1; k < 4; ++k, l = (k + 1) % 4) {
    if (!nb[k] || !nb[l]) {
      continue;
    }
    s_est = (s + (S(i + di[k], j + dj[k]) + S(i + di[l], j + dj[l]))/2)/2;
    T_new = rhr_adj(nb[k]->get_value(), nb[l]->get_value(), s_est, h, &lam);
    if (T_new < T) {
      T = T_new;
      n->set_parent_nodes({nb[k], nb[l]});
      n->set_lambda(lam);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
