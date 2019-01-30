#pragma once

#include "marcher_3d.hpp"
#include "updates.line.hpp"
#include "updates.tetra.hpp"
#include "updates.tri.hpp"

// Group dependencies:
//
// 6: IVa
// 18: I, II, III, IVb
// 26: V, VIa, VIb
//
// it may not actually be necessary to use all of these in all
// cases... This is still being investigated somewhat.

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

template <cost_func F, class base_olim3d, int num_nb>
struct abstract_olim3d:
  public marcher_3d<abstract_olim3d<F, base_olim3d, num_nb>, num_nb>
{
  static constexpr cost_func F_ = F;

  static_assert(num_nb == 6 || num_nb == 18 || num_nb == 26,
                "Number of neighbors must be 6, 18, or 26");

  using marcher_3d_t = marcher_3d<
    abstract_olim3d<F, base_olim3d, num_nb>, num_nb>;

  abstract_olim3d() { init(); }

  abstract_olim3d(vec3<int> dims, double h, no_slow_t const &):
    marcher_3d_t {dims, h, no_slow_t {}}
#if COLLECT_STATS
    , _node_stats {new updates::stats<3>[dims.product()]}
#endif
  { init(); }

  abstract_olim3d(vec3<int> dims, double h = 1,
                  std::function<double(vec3<double>)> s = static_cast<slow<3>>(s0<3>),
                  vec3<double> origin = vec3<double>::zero()):
      marcher_3d_t {dims, h, s, origin}
#if COLLECT_STATS
      , _node_stats {new updates::stats<3>[dims.product()]}
#endif
  { init(); }

  abstract_olim3d(vec3<int> dims, double h, double const * s_cache):
      marcher_3d_t {dims, h, s_cache}
#if COLLECT_STATS
      , _node_stats {new updates::stats<3>[dims.product()]}
#endif
  { init(); }

  void init();
  virtual void update_impl(int lin_hat, int * nb, int parent, double & U);

  double s_hat, s[num_nb];

#if COLLECT_STATS
  virtual ~abstract_olim3d() { delete[] _node_stats; }
  void dump_stats() const;
  void write_stats_bin(const char * path) const;
  updates::stats<3> * get_stats(int i, int j, int k) const {
    return &_node_stats[this->to_linear_index(i, j, k)];
  }
  updates::stats<3> * _stats {nullptr};
  updates::stats<3> * _node_stats {nullptr};
#endif
};

template <cost_func F, class groups>
struct olim3d_bv:
  public abstract_olim3d<F, olim3d_bv<F, groups>, groups::num_nb>
{
  static constexpr int num_nb = groups::num_nb;

  using abstract_olim3d<F, olim3d_bv<F, groups>, num_nb>::abstract_olim3d;

  void init_crtp() {}

  int lin_hat;
  int * nb;
  int parent, octant;
  int const * inds;
  bool tri_skip_list[42*8]; // TODO: compress

  void update_crtp(double & U);

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
  inline void line(int i, double & u) {
    if (this->nb[i] != -1) {
      auto u_hat = updates::line_bv<F, d>()(
        this->_U[this->nb[i]], this->s_hat, this->s[i], this->get_h());
      u = std::min(u, u_hat);
#if COLLECT_STATS
      ++this->_stats->count[0];
#endif
    }
  }

  template <int a, int b, int p0, int p1>
  inline void tri(double & u) {
    if (skip_tri<a, b>()) {
      return;
    }
    int l0 = inds[a], l1 = inds[b];
    if ((l0 == parent || l1 == parent) && nb[l0] != -1 && nb[l1] != -1) {
      auto info = updates::tri_bv<F, 3, p0, p1>()(
        this->_U[this->nb[l0]],
        this->_U[this->nb[l1]],
        this->s_hat,
        this->s[l0],
        this->s[l1],
        this->get_h());
      u = std::min(u, info.value);
#if COLLECT_STATS
      ++this->_stats->count[1];
#endif
      skip_tri<a, b>() = 1;
    }
  }

  template <int a, int b>
  inline void tri_fac(double & u) {
    if (skip_tri<a, b>()) {
      return;
    }
    int l0 = inds[a], l1 = inds[b];
    if ((l0 == parent || l1 == parent) && nb[l0] != -1 && nb[l1] != -1) {
      auto fc = this->_lin2fac[lin_hat];
      vec3<double> p_fac = fc->coords - this->to_vector_index(lin_hat);
      auto info = updates::tri<F, 3>()(
        get_p<3>(l0),
        get_p<3>(l1),
        this->_U[this->nb[l0]],
        this->_U[this->nb[l1]],
        this->s_hat,
        this->s[l0],
        this->s[l1],
        this->get_h(),
        p_fac,
        fc->s);
      u = std::min(u, info.value);
#if COLLECT_STATS
      ++this->_stats->count[1];
#endif
      skip_tri<a, b>() = 1;
    }
  }

  template <int a, int b, int c, int p0, int p1, int p2>
  inline void tetra(double & u) {
    int l0 = inds[a], l1 = inds[b], l2 = inds[c];
    if ((l0 == parent || l1 == parent || l2 == parent) &&
        this->nb[l0] != -1 && this->nb[l1] != -1 && this->nb[l2] != -1) {
      updates::info<2> info;
      double u0 = this->_U[this->nb[l0]], u1 = this->_U[this->nb[l1]],
        u2 = this->_U[this->nb[l2]], s = this->s_hat, s0 = this->s[l0],
        s1 = this->s[l1], s2 = this->s[l2], h = this->get_h();
      F_wkspc<F, 2> w;
      set_args<F>(w, u0, u1, u2, s, s0, s1, s2, h);
      cost_functor_bv<F, 3, p0, p1, p2> func {w};
      updates::tetra_bv<F, 3, p0, p1, p2>()(func, info);
      bool inbounds = info.inbounds();
      if (F == MP0 && inbounds) {
        func.set_lambda(info.lambda);
        eval_mp1_fix(func.w, s, s0, s1, s2, h, info.lambda, info.value);
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
#if COLLECT_STATS
      ++this->_stats->count[2];
#endif
    }
  }

  // TODO: this is a mess... clean it up
  template <int a, int b, int c>
  inline void tetra_fac(double & u) {
    int l0 = inds[a], l1 = inds[b], l2 = inds[c];
    if ((l0 == parent || l1 == parent || l2 == parent) &&
        this->nb[l0] != -1 && this->nb[l1] != -1 && this->nb[l2] != -1) {
      auto fc = this->_lin2fac[lin_hat];
      vec3<double> p_fac = fc->coords - this->to_vector_index(lin_hat);
      geom_fac_wkspc<2> g;
      g.init<3>(get_p<3>(l0), get_p<3>(l1), get_p<3>(l2), p_fac);
      double u0 = this->_U[this->nb[l0]], u1 = this->_U[this->nb[l1]],
        u2 = this->_U[this->nb[l2]], s = this->s_hat, s0 = this->s[l0],
        s1 = this->s[l1], s2 = this->s[l2], h = this->get_h(), s_fac = fc->s;
      F_fac_wkspc<F, 2> w;
      set_args<F>(w, g, u0, u1, u2, s, s0, s1, s2, h, s_fac);
      cost_functor_fac<F, 3, 2> func {w, g};
      updates::info<2> info;
      updates::tetra<F, 3>()(func, info);
      if (F == MP0) {
        eval_mp1_fix(func.w, s, s0, s1, s2, h, info.lambda, info.value);
      }
      u = std::min(u, info.value);
#if COLLECT_STATS
      ++this->_stats->count[2];
#endif
      skip_tri<a, b>() = 1;
      skip_tri<b, c>() = 1;
      skip_tri<a, c>() = 1;
    }
  }
};

template <class groups> using olim3d_mp0 = olim3d_bv<MP0, groups>;
template <class groups> using olim3d_mp1 = olim3d_bv<MP1, groups>;
template <class groups> using olim3d_rhr = olim3d_bv<RHR, groups>;

using olim6_groups = groups_t<0, 0, 0, 1, 0, 0, 0, 0>;
using olim6_mp0 = olim3d_mp0<olim6_groups>;
using olim6_mp1 = olim3d_mp1<olim6_groups>;
using olim6_rhr = olim3d_rhr<olim6_groups>;

using olim18_groups = groups_t<1, 0, 0, 1, 1, 0, 0, 0>;
using olim18_mp0 = olim3d_mp0<olim18_groups>;
using olim18_mp1 = olim3d_mp1<olim18_groups>;
using olim18_rhr = olim3d_rhr<olim18_groups>;

using olim26_groups = groups_t<0, 0, 0, 0, 0, 1, 0, 0>;
using olim26_mp0 = olim3d_mp0<olim26_groups>;
using olim26_mp1 = olim3d_mp1<olim26_groups>;
using olim26_rhr = olim3d_rhr<olim26_groups>;

enum LP_NORM {L1, L2, MAX};

template <cost_func F, int lp_norm, int d1, int d2>
struct olim3d_hu: public abstract_olim3d<F, olim3d_hu<F, lp_norm, d1, d2>, 26>
{
  static_assert(lp_norm == L1 || lp_norm == L2 || lp_norm == MAX,
                "Bad choice of lp norm: must be L1, L2, or MAX");
  static_assert(1 <= d1 && d1 <= 3, "d1 must satisfy 1 <= d1 <= 3");
  static_assert(1 <= d2 && d2 <= 3, "d2 must satisfy 1 <= d2 <= 3");

  using abstract_olim3d<F, olim3d_hu<F, lp_norm, d1, d2>, 26>::abstract_olim3d;

  ~olim3d_hu() {
    delete[] valid_d1;
    delete[] valid_d2;
    delete[] coplanar;
    delete[] geom_wkspcs;
    delete[] qr_wkspcs;
  }

  void init_crtp();

  int lin_hat;
  int * nb;
  int parent;
  bool * valid_d1, * valid_d2, * coplanar;
  geom_wkspc<2> * geom_wkspcs;
  qr_wkspc<3, 2> * qr_wkspcs;
  vec3<double> p0, p1, p2, p_fac;

  inline double get_p_norm(int l) const {
    if (l < 6) return 1;
    else if (l < 18) return sqrt2;
    else return sqrt3;
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

  void update_crtp(double & U);
};

using olim3d_hu_rhr = olim3d_hu<RHR, L1, 1, 2>;
using olim3d_hu_mp0 = olim3d_hu<MP0, L1, 1, 2>;
using olim3d_hu_mp1 = olim3d_hu<MP1, L1, 1, 2>;

#include "olim3d.impl.hpp"
