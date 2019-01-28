#pragma once

#include <src/config.hpp>

#include "common.hpp"
#include "vec.hpp"

template <cost_func F, bool do_adj, bool do_diag>
void
olim<F, do_adj, do_diag>::update_impl(int lin_hat, double & U)
{
  int i_hat = this->get_i(lin_hat);
  int j_hat = this->get_j(lin_hat);

  for (int k = 0; k < num_nb; ++k) {
    if (this->nb[k] != -1) {
      s[k] = this->get_s(i_hat + di<2>[k], j_hat + dj<2>[k]);
    }
  }

  if (this->is_factored(lin_hat)) {
    // TODO: this is a rough draft quality implementation of additive
    // local factoring... this can definitely be optimized

    auto fc = this->_lin2fac[lin_hat];
    double sf = fc->s;
    vec2<double> pf = {fc->i - i_hat, fc->j - j_hat};

    for (int a = 0, b = 1; a < 4; b = (++a + 1) % 4) {
      line<1>(a, U);
      tri_fac(a, b, pf, sf, U);
    }
    if (do_diag) {
      for (int a = 4, b = 0, c = 1; a < 8; ++a, c = (++b + 1) % 4) {
        line<2>(a, U);
        tri_fac(a, b, pf, sf, U);
        tri_fac(a, c, pf, sf, U);
      }
    }
  }
  else {
    for (int a = 0, b = 1; a < 4; b = (++a + 1) % 4) {
      line<1>(a, U);
      tri<0b01, 0b10>(a, b, U);
    }
    if (do_diag) {
      for (int a = 4, b = 0, c = 1; a < 8; ++a, c = (++b + 1) % 4) {
        line<2>(a, U);
        tri<0b11, 0b01>(a, b, U);
        tri<0b11, 0b10>(a, c, U);
      }
    }
  }
}
