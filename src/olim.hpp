#ifndef __OLIM_HPP__
#define __OLIM_HPP__

#include <algorithm>
#include <type_traits>

#include "marcher.hpp"
#include "node.hpp"
#include "updates.line.hpp"
#include "updates.tri.hpp"

template <cost_func F, class node, bool do_adj, bool do_diag>
struct olim:
  public marcher<
    olim<F, node, do_adj, do_diag>,
    node,
    do_diag ? 8 : 4>
{
  static constexpr int num_neighbors = do_diag ? 8 : 4;
  
  using marcher<olim<F, node, do_adj, do_diag>, node, num_neighbors>::marcher;
  static_assert(do_adj || do_diag, "error");

EIKONAL_PRIVATE:
  virtual void update_impl(node * n, double & T);

  double s_hat, s[num_neighbors];
  node * nb[num_neighbors];

  template <int d>
  inline void line(int i, double & u) {
    if (this->nb[i]) {
      u = std::min(u, updates::line_bv<F, d>()(
        this->nb[i]->get_value(), this->s_hat, this->s[i], this->get_h()));
    }
  }

  template <char p0, char p1>
  inline void tri(int i, int j, double & u) {
    if (this->nb[i] && this->nb[j]) {
      u = std::min(u, updates::tri_bv<F, 2, p0, p1>()(
        this->nb[i]->get_value(), // u0
        this->nb[j]->get_value(), // u1
        this->s_hat,              // s
        this->s[i],               // s0
        this->s[j],               // s1
        this->get_h()).value);    // h
    }
  }

  inline void tri_fac(int i, int j, double const * pf, double sf, double & u) {
    if (this->nb[i] && this->nb[j]) {
      double p0[2] = {(double) di<2>[i], (double) dj<2>[i]};
      double p1[2] = {(double) di<2>[j], (double) dj<2>[j]};
      u = std::min(u, updates::tri<F, 2>()(
        p0, p1, this->nb[i]->get_value(), this->nb[j]->get_value(),
        this->s_hat, this->s[i], this->s[j], this->get_h(), pf, sf).value);
    }
  }
};

using olim4_mp0 = olim<MP0, node, true, false>;
using olim4_mp1 = olim<MP1, node, true, false>;
using olim4_rhr = olim<RHR, node, true, false>;

using olim8_mp0 = olim<MP0, node, true, true>;
using olim8_mp1 = olim<MP1, node, true, true>;
using olim8_rhr = olim<RHR, node, true, true>;

#include "olim.impl.hpp"

#endif // __OLIM_HPP__
