#ifndef __OLIM3D_HPP__
#define __OLIM3D_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"
#if COLLECT_STATS
#    include "stats.hpp"
#endif
#include "update_rules.line_updates.hpp"
#include "update_rules.tetra_updates.hpp"
#include "update_rules.tri_updates.hpp"

// dependencies:
//
// 6: IVa
// 18: I, II, III, IVb
// 26: V, VIa, VIb

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
  static constexpr int nneib = V || VI_a || VI_b ? 26 :
    (I || II || III || IV_b ? 18 : 6);
};

template <class base_olim3d, class node, class line_updates, class tri_updates,
          class tetra_updates, int num_neighbors>
struct abstract_olim3d:
  public marcher_3d<
    abstract_olim3d<
      base_olim3d, node, line_updates, tri_updates, tetra_updates,
      num_neighbors>,
    node>,
  public line_updates,
  public tri_updates,
  public tetra_updates
{
  using marcher_3d_t = marcher_3d<
    abstract_olim3d<
      base_olim3d, node, line_updates, tri_updates, tetra_updates,
      num_neighbors>,
    node>;

  static constexpr int nneib = nneib;

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
  virtual void visit_neighbors_impl(abstract_node * n);
  virtual void update_impl(int i, int j, int k, int parent, double & T);
  void init();
#if COLLECT_STATS
  olim3d_node_stats * _node_stats {nullptr};
#endif
};

template <
  class node, class line_updates, class tri_updates, class tetra_updates,
  class groups>
struct olim3d_bv:
  public abstract_olim3d<
    olim3d_bv<node, line_updates, tri_updates, tetra_updates, groups>,
    node, line_updates, tri_updates, tetra_updates, groups::nneib>,
  public groups
{
  using abstract_olim3d<
    olim3d_bv<node, line_updates, tri_updates, tetra_updates, groups>,
    node, line_updates, tri_updates, tetra_updates,
    groups::nneib>::abstract_olim3d;

  void update_crtp(int i, int j, int k, int parent, double & T);
};

template <class groups>
using olim3d_mp0 = olim3d_bv<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates,
  update_rules::mp0_tetra_updates_bv,
  groups>;

template <class groups>
using olim3d_mp1 = olim3d_bv<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates,
  update_rules::mp1_tetra_updates_bv,
  groups>;

template <class groups>
using olim3d_rhr = olim3d_bv<
  node_3d,
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates,
  update_rules::rhr_tetra_updates_bv,
  groups>;

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

template <
  class node, class line_updates, class tri_updates, class tetra_updates,
  int lp_norm, int d1, int d2>
struct olim3d_hu:
  public abstract_olim3d<
    olim3d_hu<node, line_updates, tri_updates, tetra_updates, lp_norm, d1, d2>,
    node, line_updates, tri_updates, tetra_updates, 26>
{
  static_assert(lp_norm == L1 || lp_norm == L2 || lp_norm == MAX,
                "Bad choice of lp norm: must be L1, L2, or MAX");
  static_assert(1 <= d1 && d1 <= 3, "d1 must satisfy 1 <= d1 <= 3");
  static_assert(1 <= d2 && d2 <= 3, "d2 must satisfy 1 <= d2 <= 3");

  using abstract_olim3d<
    olim3d_hu<node, line_updates, tri_updates, tetra_updates, lp_norm, d1, d2>,
    node, line_updates, tri_updates, tetra_updates, 26>::abstract_olim3d;

  void update_crtp(int i, int j, int k, int parent, double & T);
};

using olim3d_hu_rhr = olim3d_hu<
  node_3d,
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates,
  update_rules::rhr_tetra_updates,
  L1, 1, 2>;

using olim3d_hu_mp0 = olim3d_hu<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates,
  update_rules::mp0_tetra_updates,
  L1, 1, 2>;

using olim3d_hu_mp1 = olim3d_hu<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates,
  update_rules::mp1_tetra_updates,
  L1, 1, 2>;

#include "olim3d.impl.hpp"

#endif // __OLIM3D_HPP__
