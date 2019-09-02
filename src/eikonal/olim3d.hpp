#pragma once

#include "line.hpp"
#include "marcher.hpp"
#include "tetra.hpp"
#include "tri.hpp"

// Group dependencies:
//
// 6: IVa
// 18: I, II, III, IVb
// 26: V, VIa, VIb
//
// it may not actually be necessary to use all of these in all
// cases... This is still being investigated somewhat.

namespace eikonal {

template <bool I, bool II, bool III, bool IV_a, bool IV_b, bool V,
          bool VI_a, bool VI_b>
struct groups_t {
  static constexpr bool group_I = I;
  static constexpr bool group_II = II;
  static constexpr bool group_III = III;
  static constexpr bool group_IV_a = IV_a;
  static constexpr bool group_IV_b = IV_b;
  static constexpr bool group_V = V;
  static constexpr bool group_VI_a = VI_a;
  static constexpr bool group_VI_b = VI_b;
  static constexpr bool do_line2_updates = I || II || III || IV_b || V || VI_b;
  static constexpr bool do_line3_updates = V || VI_a || VI_b;
  static constexpr bool do_tri11_updates = II || III || IV_a || VI_a;
  static constexpr bool do_tri12_updates = I || II || III || V;
  static constexpr bool do_tri13_updates = V;
  static constexpr bool do_tri22_updates = I || II || III || IV_b || VI_b;
  static constexpr bool do_tri23_updates = V || VI_b;
  static constexpr int num_nb = V || VI_a || VI_b ? 26 :
    (I || II || III || IV_b ? 18 : 6);
};

template <cost_func F, class groups, ordering ord>
struct olim3d_bv:
  public marcher<olim3d_bv<F, groups, ord>, 3, groups::num_nb, ord>
{
  static constexpr cost_func F_ = F;

  using marcher<olim3d_bv<F, groups, ord>, 3, groups::num_nb, ord>::marcher;

  int octant;
  int const * inds;
  bool tri_skip_list[42*8]; // TODO: compress

  void update_impl(int lin_hat, int const * nb, int parent, double & U);

  // TODO: see comment in olim.hpp
  inline double line(double l0, double u0, double s, double s0, double h) const {
    return eikonal::line<F>()(l0, u0, s, s0, h);
  }

OLIM_PRIVATE:

  inline void reset_tri_skip_list() {
    memset(tri_skip_list, 0, sizeof(tri_skip_list));
  }

  template <int i, int j>
  inline bool & skip_tri() {
    int m = i < j ? i : j, M = i < j ? j : i;
    return tri_skip_list[42*octant + 7*m + M];
  }

  template <int d>
  inline void line(int lin_hat, int const * nb, int i, double & U) {
    if (nb[i] != -1) {
      U = fmin(U, line_bv<F, d>()(
        this->_U[nb[i]],
        this->_s[lin_hat],
        this->_s[nb[i]],
        this->get_h()));
    }
  }

  template <int a, int b, int p0, int p1>
  inline void tri(int lin_hat, int const * nb, int parent, double & u) {
    if (skip_tri<a, b>()) {
      return;
    }
    int l0 = inds[a], l1 = inds[b];
    if ((l0 == parent || l1 == parent) && nb[l0] != -1 && nb[l1] != -1) {
      auto info = tri_bv<F, 3, p0, p1>()(
        this->_U[nb[l0]],
        this->_U[nb[l1]],
        this->_s[lin_hat],
        this->_s[nb[l0]],
        this->_s[nb[l1]],
        this->get_h());
      u = std::min(u, info.value);
      skip_tri<a, b>() = 1;
    }
  }

  template <int a, int b>
  inline void tri_fac(int lin_hat, int const * nb, int parent, double & u) {
    if (skip_tri<a, b>()) {
      return;
    }
    int l0 = inds[a], l1 = inds[b];
    if ((l0 == parent || l1 == parent) && nb[l0] != -1 && nb[l1] != -1) {
      auto src = this->_fac_srcs[lin_hat];
      auto info = eikonal::tri<F, 3>()(
        get_p<3>(l0),
        get_p<3>(l1),
        this->_U[nb[l0]],
        this->_U[nb[l1]],
        this->_s[lin_hat],
        this->_s[nb[l0]],
        this->_s[nb[l1]],
        this->get_h(),
        vec3<double> {src->x} - this->to_vector_index(lin_hat),
        src->s);
      u = std::min(u, info.value);
      skip_tri<a, b>() = 1;
    }
  }

  template <int a, int b, int c, int p0, int p1, int p2>
  inline void tetra(int lin_hat, int const * nb, int parent, double & u) {
    int l0 = inds[a], l1 = inds[b], l2 = inds[c];
    if ((l0 == parent || l1 == parent || l2 == parent) &&
        nb[l0] != -1 && nb[l1] != -1 && nb[l2] != -1) {
      update_info<2> info;
      F_wkspc<F, 2> w;
      set_args<F>(
        w,
        this->_U[nb[l0]],
        this->_U[nb[l1]],
        this->_U[nb[l2]],
        this->_s[lin_hat],
        this->_s[nb[l0]],
        this->_s[nb[l1]],
        this->_s[nb[l2]],
        this->get_h());
      cost_functor_bv<F, 3, p0, p1, p2> func {w};
      tetra_bv<F, 3, p0, p1, p2>()(func, info);
      bool inbounds = info.inbounds();
      if (F == MP0 && inbounds) {
        func.set_lambda(info.lambda);
        eval_mp1_fix(
          func.w,
          this->_s[lin_hat],
          this->_s[nb[l0]],
          this->_s[nb[l1]],
          this->_s[nb[l2]],
          this->get_h(),
          info.lambda,
          info.value);
      }
      u = std::min(u, info.value);
      if (F == MP1 || inbounds) {
        skip_tri<a, b>() = 1;
        skip_tri<b, c>() = 1;
        skip_tri<a, c>() = 1;
      } else if (info.finite_lambda()) { // (F == MP0 || F == RHR) && !inbounds
        auto const & lam = info.lambda;
        if (lam[0] < 0 && lam[1] < 0) {
          skip_tri<b, c>() = 1;
        } else if (lam[0] + lam[1] > 1) {
          if (lam[0] > 0 && lam[1] > 0) {
            skip_tri<a, b>() = 1;
            skip_tri<a, c>() = 1;
          } else if (lam[0] < 0) {
            skip_tri<a, c>() = 1;
          } else if (lam[1] < 0) {
            skip_tri<a, b>() = 1;
          }
        } else if (lam[0] < 0) {
          skip_tri<a, b>() = 1;
          skip_tri<b, c>() = 1;
        } else if (lam[1] < 0) {
          skip_tri<b, c>() = 1;
          skip_tri<a, c>() = 1;
        }
      }
    }
  }

  // TODO: this is a mess... clean it up
  template <int a, int b, int c>
  inline void tetra_fac(int lin_hat, int const * nb, int parent, double & u) {
    int l0 = inds[a], l1 = inds[b], l2 = inds[c];
    if ((l0 == parent || l1 == parent || l2 == parent) &&
        nb[l0] != -1 && nb[l1] != -1 && nb[l2] != -1) {
      auto src = this->_fac_srcs[lin_hat];
      auto p_fac = vec3<double> {src->x} - this->to_vector_index(lin_hat);
      geom_fac_wkspc<2> g;
      g.init<3>(get_p<3>(l0), get_p<3>(l1), get_p<3>(l2), p_fac);
      F_fac_wkspc<F, 2> w;
      set_args<F>(
        w,
        g,
        this->_U[nb[l0]],
        this->_U[nb[l1]],
        this->_U[nb[l2]],
        this->_s[lin_hat],
        this->_s[nb[l0]],
        this->_s[nb[l1]],
        this->_s[nb[l2]],
        this->get_h(),
        src->s);
      cost_functor_fac<F, 3, 2> func {w, g};
      update_info<2> info;
      eikonal::tetra<F, 3>()(func, info);
      if (F == MP0) {
        eval_mp1_fix(
          func.w,
          this->_s[lin_hat],
          this->_s[nb[l0]],
          this->_s[nb[l1]],
          this->_s[nb[l2]],
          this->get_h(),
          info.lambda,
          info.value);
      }
      u = std::min(u, info.value);
      skip_tri<a, b>() = 1;
      skip_tri<b, c>() = 1;
      skip_tri<a, c>() = 1;
    }
  }
};

using olim6_groups = groups_t<0, 0, 0, 1, 0, 0, 0, 0>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim6_mp0 = olim3d_bv<MP0, olim6_groups, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim6_mp1 = olim3d_bv<MP1, olim6_groups, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim6_rhr = olim3d_bv<RHR, olim6_groups, ord>;

using olim18_groups = groups_t<1, 0, 0, 1, 1, 0, 0, 0>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim18_mp0 = olim3d_bv<MP0, olim18_groups, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim18_mp1 = olim3d_bv<MP1, olim18_groups, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim18_rhr = olim3d_bv<RHR, olim18_groups, ord>;

using olim26_groups = groups_t<0, 0, 0, 0, 0, 1, 0, 0>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim26_mp0 = olim3d_bv<MP0, olim26_groups, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim26_mp1 = olim3d_bv<MP1, olim26_groups, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim26_rhr = olim3d_bv<RHR, olim26_groups, ord>;

enum LP_NORM {L1, L2, MAX};

template <cost_func F, int lp_norm, int d1, int d2, ordering ord>
struct olim3d_hu:
  public marcher<olim3d_hu<F, lp_norm, d1, d2, ord>, 3, 26, ord>
{
  static_assert(lp_norm == L1 || lp_norm == L2 || lp_norm == MAX,
                "Bad choice of lp norm: must be L1, L2, or MAX");
  static_assert(1 <= d1 && d1 <= 3, "d1 must satisfy 1 <= d1 <= 3");
  static_assert(1 <= d2 && d2 <= 3, "d2 must satisfy 1 <= d2 <= 3");

  static constexpr cost_func F_ = F;

  using ivec = vec<int, 3>;
  using fvec = vec<double, 3>;

  using marcher_t = marcher<olim3d_hu<F, lp_norm, d1, d2, ord>, 3, 26, ord>;

  // TODO: define an init_crtp that can be called from marcher and get
  // rid of all these explicit calls to marcher's constructors

  olim3d_hu(ivec dims, double h): marcher_t {dims, h} {
    init();
  }

  olim3d_hu(ivec dims, double h, no_slow_t const &):
    marcher_t {dims, h, no_slow_t {}}
  {
    init();
  }

  olim3d_hu(ivec dims, double h, double const * s): marcher_t {dims, h, s} {
    init();
  }

  olim3d_hu(ivec dims, double h, slow<3> s, fvec origin = fvec::zero()):
    marcher_t {dims, h, s, origin}
  {
    init();
  }

  ~olim3d_hu() {
    delete[] valid_d1;
    delete[] valid_d2;
    delete[] coplanar;
    delete[] geom_wkspcs;
    delete[] qr_wkspcs;
  }

  void init();

  bool * valid_d1, * valid_d2, * coplanar;
  geom_wkspc<2> * geom_wkspcs;
  qr_wkspc<3, 2> * qr_wkspcs;
  vec3<double> p0, p1, p2, p_fac;

  inline double get_p_norm(int l) const {
    if (l < 6) return 1;
    else if (l < 18) return int_sqrt(2);
    else return int_sqrt(3);
  };

  inline int to_nb_linear_index(vec2<int> inds) {
    return 26*inds[0] + inds[1];
  }

  inline int to_nb_linear_index(vec3<int> inds) {
    return 26*(26*inds[0] + inds[1]) + inds[2];
  }

  inline geom_wkspc<2> & get_geom_wkspc(vec3<int> inds) {
    return geom_wkspcs[to_nb_linear_index(inds)];
  }

  inline qr_wkspc<3, 2> & get_qr_wkspc(vec3<int> inds) {
    return qr_wkspcs[to_nb_linear_index(inds)];
  }

  inline bool & is_valid_d1(vec2<int> inds) {
    return valid_d1[to_nb_linear_index(inds)];
  }

  inline bool & is_valid_d2(vec2<int> inds) {
    return valid_d2[to_nb_linear_index(inds)];
  }

  inline bool & is_coplanar(vec3<int> inds) {
    return coplanar[to_nb_linear_index(inds)];
  }

  // TODO: see comment in olim.hpp
  inline double line(double l0, double u0, double s, double s0, double h) {
    return eikonal::line<F>()(l0, u0, s, s0, h);
  }

  void update_impl(int lin_hat, int const * nb, int parent, double & U);
};

template <ordering ord = ordering::COLUMN_MAJOR>
using olim3d_rhr = olim3d_hu<RHR, L1, 1, 2, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim3d_mp0 = olim3d_hu<MP0, L1, 1, 2, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using olim3d_mp1 = olim3d_hu<MP1, L1, 1, 2, ord>;

}

#include "olim3d.impl.hpp"
