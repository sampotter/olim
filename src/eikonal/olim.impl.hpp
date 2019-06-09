#pragma once

#include <src/config.hpp>

#include "../common.hpp"
#include "../vec.hpp"

template <cost_func F, bool do_adj, bool do_diag, ordering ord>
void
olim<F, do_adj, do_diag, ord>::update_impl(
  int lin_hat, int const * nb, int parent, double & U)
{
  (void) parent;

  vec2<int> const inds = this->to_vector_index(lin_hat);

  if (this->is_factored(lin_hat)) {
    // TODO: this is a rough draft quality implementation of additive
    // local factoring... this can definitely be optimized

    auto fc = this->_lin2fac[lin_hat];
    double sf = fc->s;
    vec2<double> pf = fc->coords - inds;

    for (int a = 0, b = 1; a < 4; b = (++a + 1) % 4) {
      line<1>(lin_hat, nb, a, U);
      tri_fac(lin_hat, nb, a, b, pf, sf, U);
    }
    if (do_diag) {
      for (int a = 4, b = 0, c = 1; a < 8; ++a, c = (++b + 1) % 4) {
        line<2>(lin_hat, nb, a, U);
        tri_fac(lin_hat, nb, a, b, pf, sf, U);
        tri_fac(lin_hat, nb, a, c, pf, sf, U);
      }
    }
  }
  else {
    for (int a = 0, b = 1; a < 4; b = (++a + 1) % 4) {
      line<1>(lin_hat, nb, a, U);
      tri<0b01, 0b10>(lin_hat, nb, a, b, U);
    }
    if (do_diag) {
      for (int a = 4, b = 0, c = 1; a < 8; ++a, c = (++b + 1) % 4) {
        line<2>(lin_hat, nb, a, U);
        tri<0b11, 0b01>(lin_hat, nb, a, b, U);
        tri<0b11, 0b10>(lin_hat, nb, a, c, U);
      }
    }
  }
}
