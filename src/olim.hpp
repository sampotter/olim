#pragma once

#include <type_traits>

#include "marcher.hpp"
#include "updates.line.hpp"
#include "updates.tri.hpp"

template <cost_func F, bool do_adj, bool do_diag>
struct olim: public marcher<olim<F, do_adj, do_diag>, 2, do_diag ? 8 : 4>
{
  static constexpr cost_func F_ = F;
  static constexpr int num_nb = do_diag ? 8 : 4;

  using marcher<olim<F, do_adj, do_diag>, 2, num_nb>::marcher;
  static_assert(do_adj || do_diag, "error");

OLIM_PRIVATE:
  virtual void update_impl(int lin_hat, int const * nb, int parent, double & U);

  template <int d>
  inline void line(int lin_hat, int const * nb, int i, double & u) {
    if (nb[i] != -1) {
      u = std::min(u, updates::line_bv<F, d>()(
        this->_U[nb[i]],
        this->_s[lin_hat],
        this->_s[nb[i]],
        this->get_h()));
    }
  }

  template <int p0, int p1>
  inline void tri(int lin_hat, int const * nb, int i, int j, double & u) {
    if (nb[i] != -1 && nb[j] != -1) {
      u = std::min(u, updates::tri_bv<F, 2, p0, p1>()(
        this->_U[nb[i]],
        this->_U[nb[j]],
        this->_s[lin_hat],
        this->_s[nb[i]],
        this->_s[nb[j]],
        this->get_h()
      ).value);
    }
  }

  inline void tri_fac(int lin_hat, int const * nb,
                      int i, int j, vec2<double> const & pf, double sf,
                      double & u) {
    if (nb[i] != -1 && nb[j] != -1) {
      vec2<double> p0 = get_p<2>(i), p1 = get_p<2>(j);
      u = std::min(u, updates::tri<F, 2>()(
        p0,
        p1,
        this->_U[nb[i]],
        this->_U[nb[j]],
        this->_s[lin_hat],
        this->_s[i],
        this->_s[j],
        this->get_h(),
        pf,
        sf
      ).value);
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
