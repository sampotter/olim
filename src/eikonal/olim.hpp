#pragma once

// TODO: remove this
#include <type_traits>

#include "line.hpp"
#include "marcher.hpp"
#include "tri.hpp"

template <cost_func F, bool do_adj, bool do_diag, ordering ord>
struct olim: public marcher<
  olim<F, do_adj, do_diag, ord>, 2, do_diag ? 8 : 4, ord>
{
  static constexpr cost_func F_ = F;
  static constexpr int num_nb = do_diag ? 8 : 4;

  using marcher<olim<F, do_adj, do_diag, ord>, 2, num_nb, ord>::marcher;
  static_assert(do_adj || do_diag, "error");

OLIM_PRIVATE:
  virtual void update_impl(int lin_hat, int const * nb, int parent, double & U);

  // TODO: this is used when we add non-grid-aligned points... there
  // is probably a nicer, more flexible way of doing this, but it
  // works for now.
  inline double line(double l0, double u0, double s, double s0, double h) const {
    return eikonal::line<F>()(l0, u0, s, s0, h);
  }

  template <int d>
  inline void line(int lin_hat, int const * nb, int i, double & u) {
    if (nb[i] != -1) {
      u = std::min(u, eikonal::line_bv<F, d>()(
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

template <ordering ord = ordering::COLUMN_MAJOR>
using olim4_mp0 = olim<MP0, true, false, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim4_mp1 = olim<MP1, true, false, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim4_rhr = olim<RHR, true, false, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim8_mp0 = olim<MP0, false, true, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim8_mp1 = olim<MP1, false, true, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim8_rhr = olim<RHR, false, true, ord>;

#include "olim.impl.hpp"
