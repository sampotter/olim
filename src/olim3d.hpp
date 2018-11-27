#ifndef __OLIM3D_HPP__
#define __OLIM3D_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"
#if COLLECT_STATS
#    include "stats.hpp"
#endif
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
  static constexpr int num_neighbors = V || VI_a || VI_b ? 26 :
    (I || II || III || IV_b ? 18 : 6);
};

template <class base_olim3d, class node, int num_neighbors>
struct abstract_olim3d:
  public marcher_3d<
    abstract_olim3d<base_olim3d, node, num_neighbors>,
    node,
    num_neighbors>
{
  static constexpr int num_nb = num_neighbors;
  static_assert(num_nb == 6 || num_nb == 18 || num_nb == 26,
                "Number of neighbors must be 6, 18, or 26");

  using marcher_3d_t = marcher_3d<
    abstract_olim3d<base_olim3d, node, num_neighbors>,
    node,
    num_neighbors>;

  abstract_olim3d() {}

  abstract_olim3d(int height, int width, int depth, double h,
                  no_speed_func_t const &):
      marcher_3d_t {height, width, depth, h, no_speed_func_t {}}
#if COLLECT_STATS
    , _node_stats {new olim3d_node_stats[height*width*depth]} {}
#else
    {}
#endif

  abstract_olim3d(int height, int width, int depth, double h = 1,
                  std::function<double(double, double, double)> speed =
                    static_cast<speed_func_3d>(default_speed_func),
                  double x0 = 0.0, double y0 = 0.0, double z0 = 0.0):
      marcher_3d_t {height, width, depth, h, speed, x0, y0, z0}
#if COLLECT_STATS
    , _node_stats {new olim3d_node_stats[height*width*depth]} {}
#else
    {}
#endif

  abstract_olim3d(int height, int width, int depth, double h,
                  double const * s_cache):
      marcher_3d_t {height, width, depth, h, s_cache}
#if COLLECT_STATS
    , _node_stats {new olim3d_node_stats[height*width*depth]} {}
#else
    {}
#endif

#if COLLECT_STATS
  virtual ~abstract_olim3d() { delete[] _node_stats; }
  void dump_stats() const;
  olim3d_node_stats & get_node_stats(int i, int j, int k);
  olim3d_node_stats const & get_node_stats(int i, int j, int k) const;
#endif

EIKONAL_PRIVATE:
  virtual void update_impl(node * n, node ** nb, int parent, double & T);
  void init();

  double s_hat, s[num_neighbors];

#if COLLECT_STATS
  olim3d_node_stats * _node_stats {nullptr};
#endif
};

template <cost_func F, class node, class groups>
struct olim3d_bv:
  public abstract_olim3d<olim3d_bv<F, node, groups>, node, groups::num_neighbors>
{
  static constexpr int num_neighbors = groups::num_neighbors;
  
  using abstract_olim3d<
    olim3d_bv<F, node, groups>, node, num_neighbors>::abstract_olim3d;

EIKONAL_PROTECTED:
  void update_crtp(double & T);

  node * n;
  node ** nb;
  int parent, octant;
  int const * inds;
  bool tri_skip_list[42*8]; // TODO: compress

  inline void reset_tri_skip_list() {
    for (int i = 0; i < 42*8; ++i) {
      tri_skip_list[i] = 0;
    }
  }

  template <int i, int j>
  inline bool & skip_tri() {
    int m = i < j ? i : j, M = i < j ? j : i;
    return tri_skip_list[42*octant + 7*m + M];
  }

  template <int d>
  inline void line(int i, double & u) {
    if (this->nb[i]) {
      auto u_hat = updates::line_bv<F, d>()(
        this->nb[i]->get_value(), this->s_hat, this->s[i], this->get_h());
#if TRACK_PARENTS
#  error Not implemented yet!
#else
      u = std::min(u, u_hat);
#endif
#if COLLECT_STATS
      node_stats.add_line_update(d);
#endif
    }
  }

  template <int a, int b, char p0, char p1>
  inline void tri(double & u) {
    if (skip_tri<a, b>()) {
      return;
    }
    int l0 = inds[a], l1 = inds[b];
    if ((l0 == parent || l1 == parent) && nb[l0] && nb[l1]) {
      auto info = updates::tri_bv<F, 3, p0, p1>()(
        this->nb[l0]->get_value(),
        this->nb[l1]->get_value(),
        this->s_hat,
        this->s[l0],
        this->s[l1],
        this->get_h());
#if TRACK_PARENTS
#  error Not implemented yet!
#else
      u = std::min(u, info.value);
#endif
#if COLLECT_STATS
      node_stats.add_tri_update(p0, p1, info);
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
    if ((l0 == parent || l1 == parent) && nb[l0] && nb[l1]) {
      auto n_fac = static_cast<node *>(n->get_fac_parent());
      int i_fac = n_fac->get_i();
      int j_fac = n_fac->get_j();
      int k_fac = n_fac->get_k();
      double p0[3] = {(double)di<3>[l0], (double)dj<3>[l0], (double)dk<3>[l0]};
      double p1[3] = {(double)di<3>[l1], (double)dj<3>[l1], (double)dk<3>[l1]};
      double p_fac[3] = {
        (double) (i_fac - n->get_i()),
        (double) (j_fac - n->get_j()),
        (double) (k_fac - n->get_k())
      };
      auto info = updates::tri<F, 3>()(
        p0,
        p1,
        this->nb[l0]->get_value(),
        this->nb[l1]->get_value(),
        this->s_hat,
        this->s[l0],
        this->s[l1],
        this->get_h(),
        p_fac,
        this->get_speed(i_fac, j_fac, k_fac));
#if TRACK_PARENTS
#  error Not implemented yet!
#else
      u = std::min(u, info.value);
#endif
#if COLLECT_STATS
#  error Not implemented yet!
#endif
      skip_tri<a, b>() = 1;
    }
  }

  template <int a, int b, int c, char p0, char p1, char p2>
  inline void tetra(double & u) {
    int l0 = inds[a], l1 = inds[b], l2 = inds[c];
    if ((l0 == parent || l1 == parent || l2 == parent) &&
        this->nb[l0] && this->nb[l1] && this->nb[l2]) {
      auto info = updates::tetra_bv<F, 3, p0, p1, p2>()(
        this->nb[l0]->get_value(),
        this->nb[l1]->get_value(),
        this->nb[l2]->get_value(),
        this->s_hat,
        this->s[l0],
        this->s[l1],
        this->s[l2],
        this->get_h());
#if TRACK_PARENTS
#  error Not implemented yet!
#else
      u = std::min(u, info.value);
#endif
#if COLLECT_STATS
      node_stats.add_tetra_update(p0, p1, p2, info);
#endif
      skip_tri<a, b>() = 1;
      skip_tri<b, c>() = 1;
      skip_tri<a, c>() = 1;
    }
  }

  template <int a, int b, int c>
  inline void tetra_fac(double & u) {
    int l0 = inds[a], l1 = inds[b], l2 = inds[c];
    if ((l0 == parent || l1 == parent || l2 == parent) &&
        this->nb[l0] && this->nb[l1] && this->nb[l2]) {
      auto n_fac = static_cast<node *>(n->get_fac_parent());
      int i_fac = n_fac->get_i();
      int j_fac = n_fac->get_j();
      int k_fac = n_fac->get_k();
      double p0[3] = {(double)di<3>[l0], (double)dj<3>[l0], (double)dk<3>[l0]};
      double p1[3] = {(double)di<3>[l1], (double)dj<3>[l1], (double)dk<3>[l1]};
      double p2[3] = {(double)di<3>[l2], (double)dj<3>[l2], (double)dk<3>[l2]};
      double p_fac[3] = {
        (double) (i_fac - n->get_i()),
        (double) (j_fac - n->get_j()),
        (double) (k_fac - n->get_k())
      };
      auto info = updates::tetra<F, 3>()(
        p0,
        p1,
        p2,
        this->nb[l0]->get_value(),
        this->nb[l1]->get_value(),
        this->nb[l2]->get_value(),
        this->s_hat,
        this->s[l0],
        this->s[l1],
        this->s[l2],
        this->get_h(),
        p_fac,
        this->get_speed(i_fac, j_fac, k_fac));
#if TRACK_PARENTS
#  error Not implemented yet!
#else
      u = std::min(u, info.value);
#endif
#if COLLECT_STATS
#  error Not implemented yet!
#endif
      skip_tri<a, b>() = 1;
      skip_tri<b, c>() = 1;
      skip_tri<a, c>() = 1;
    }
  }
};

template <class groups> using olim3d_mp0 = olim3d_bv<MP0, node_3d, groups>;
template <class groups> using olim3d_mp1 = olim3d_bv<MP1, node_3d, groups>;
template <class groups> using olim3d_rhr = olim3d_bv<RHR, node_3d, groups>;

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

template <cost_func F, class node, int lp_norm, int d1, int d2>
struct olim3d_hu:
  public abstract_olim3d<olim3d_hu<F, node, lp_norm, d1, d2>, node, 26>
{
  static_assert(lp_norm == L1 || lp_norm == L2 || lp_norm == MAX,
                "Bad choice of lp norm: must be L1, L2, or MAX");
  static_assert(1 <= d1 && d1 <= 3, "d1 must satisfy 1 <= d1 <= 3");
  static_assert(1 <= d2 && d2 <= 3, "d2 must satisfy 1 <= d2 <= 3");

  using abstract_olim3d<
    olim3d_hu<F, node, lp_norm, d1, d2>, node, 26>::abstract_olim3d;

  double s_hat, s[26];
  node * n;
  node ** nb;
  int parent, octant;
  int const * inds;
  int tri_skip_list[42*8];

  void update_crtp(double & T);
};

using olim3d_hu_rhr = olim3d_hu<RHR, node_3d, L1, 1, 2>;
using olim3d_hu_mp0 = olim3d_hu<MP0, node_3d, L1, 1, 2>;
using olim3d_hu_mp1 = olim3d_hu<MP1, node_3d, L1, 1, 2>;

#include "olim3d.impl.hpp"

#endif // __OLIM3D_HPP__
