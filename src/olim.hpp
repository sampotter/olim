#pragma once

#include <type_traits>

#include "marcher.hpp"
#include "updates.line.hpp"
#include "updates.tri.hpp"

template <cost_func F, bool do_adj, bool do_diag>
struct olim: public marcher<olim<F, do_adj, do_diag>, do_diag ? 8 : 4>
{
  static constexpr cost_func F_ = F;
  static constexpr int num_nb = do_diag ? 8 : 4;

  using marcher<olim<F, do_adj, do_diag>, num_nb>::marcher;
  static_assert(do_adj || do_diag, "error");

  double s_hat, s[num_nb];
  int nb[num_nb];

OLIM_PRIVATE:
  virtual void update_impl(int lin_hat, double & U);

  template <int d>
  inline void line(int i, double & u) {
    if (this->nb[i] != -1) {
      u = std::min(u, updates::line_bv<F, d>()(
        this->_U[this->nb[i]], this->s_hat, this->s[i], this->get_h()));
    }
  }

  template <int p0, int p1>
  inline void tri(int i, int j, double & u) {
    if (this->nb[i] != -1 && this->nb[j] != -1) {
      u = std::min(u, updates::tri_bv<F, 2, p0, p1>()(
        this->_U[this->nb[i]], // u0
        this->_U[this->nb[j]], // u1
        this->s_hat,           // s
        this->s[i],            // s0
        this->s[j],            // s1
        this->get_h()).value); // h
    }
  }

  inline void tri_fac(int i, int j, vec2<double> const & pf, double sf, double & u) {
    if (this->nb[i] != -1 && this->nb[j] != -1) {
      vec2<double> p0 = {(double) di<2>[i], (double) dj<2>[i]};
      vec2<double> p1 = {(double) di<2>[j], (double) dj<2>[j]};
      u = std::min(u, updates::tri<F, 2>()(
        p0, p1, this->_U[this->nb[i]], this->_U[this->nb[j]],
        this->s_hat, this->s[i], this->s[j], this->get_h(), pf, sf).value);
    }
  }
};

using olim4_mp0 = olim<MP0, true, false>;
using olim4_mp1 = olim<MP1, true, false>;
using olim4_rhr = olim<RHR, true, false>;

using olim8_mp0 = olim<MP0, true, true>;
using olim8_mp1 = olim<MP1, true, true>;
using olim8_rhr = olim<RHR, true, true>;

#include "olim.impl.hpp"
