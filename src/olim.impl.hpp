#ifndef __OLIM_IMPL_HPP__
#define __OLIM_IMPL_HPP__

#include <src/config.hpp>

#include "common.hpp"
#include "vecmath.hpp"

#define __di(k) di<2>[k]
#define __dj(k) dj<2>[k]

#define P01 1
#define P10 2
#define P11 3

template <cost_func F, class node, bool do_adj, bool do_diag>
void
olim<F, node, do_adj, do_diag>::update_impl(node * n, double & T)
{
  int i_hat = n->get_i(), j_hat = n->get_j();

  for (int k = 0; k < num_neighbors; ++k) {
    if (this->nb[k]) {
      s[k] = this->get_speed(i_hat + __di(k), j_hat + __dj(k));
    }
  }

  if (n->is_factored()) {
    // TODO: this is a rough draft quality implementation of additive
    // local factoring... this can definitely be optimized

    auto fc = n->get_fac_center();
    double sf = fc->s, pf[2] = {fc->i - i_hat, fc->j - j_hat};

    for (int a = 0, b = 1; a < 4; b = (++a + 1) % 4) {
      line<1>(a, T);
      tri_fac(a, b, pf, sf, T);
    }
    if (do_diag) {
      for (int a = 4, b = 0, c = 1; a < 8; ++a, c = (++b + 1) % 4) {
        line<2>(a, T);
        tri_fac(a, b, pf, sf, T);
        tri_fac(a, c, pf, sf, T);
      }
    }
  }
  else {
    for (int a = 0, b = 1; a < 4; b = (++a + 1) % 4) {
      line<1>(a, T);
      tri<P01, P10>(a, b, T);
    }
    if (do_diag) {
      for (int a = 4, b = 0, c = 1; a < 8; ++a, c = (++b + 1) % 4) {
        line<2>(a, T);
        tri<P11, P01>(a, b, T);
        tri<P11, P10>(a, c, T);
      }
    }
  }
}

#undef __di
#undef __dj

#undef P01
#undef P10
#undef P11
#undef LINE
#undef DO_LINE
#undef TRI
#undef DO_TRI
#undef TRI_FAC
#undef DO_TRI_FAC

#endif // __OLIM_IMPL_HPP__
