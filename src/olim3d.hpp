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
          class tetra_updates, int nneib>
struct abstract_olim3d: public marcher_3d<node>, public line_updates,
                        public tri_updates, public tetra_updates
{
  abstract_olim3d(int height, int width, int depth, double h,
                  no_speed_func_t const &):
      marcher_3d<node> {height, width, depth, h, no_speed_func_t {}}
#if COLLECT_STATS
    , _node_stats {new olim3d_node_stats[height*width*depth]} {}
#else
    {}
#endif

  abstract_olim3d(int height, int width, int depth, double h = 1,
                  std::function<double(double, double, double)> speed =
                    static_cast<speed_func_3d>(default_speed_func),
                  double x0 = 0.0, double y0 = 0.0, double z0 = 0.0):
      marcher_3d<node> {height, width, depth, h, speed, x0, y0, z0}
#if COLLECT_STATS
    , _node_stats {new olim3d_node_stats[height*width*depth]} {}
#else
    {}
#endif

  abstract_olim3d(int height, int width, int depth, double h,
                  double const * s_cache):
      marcher_3d<node> {height, width, depth, h, s_cache}
#if COLLECT_STATS
    , _node_stats {new olim3d_node_stats[height*width*depth]} {}
#else
    {}
#endif

#if COLLECT_STATS
  virtual ~olim3d() { delete[] _node_stats; }
  void dump_stats() const;
#endif

EIKONAL_PROTECTED:
  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb);

EIKONAL_PRIVATE:
  virtual void stage_neighbors_impl(abstract_node * n);
  virtual void update_impl(int i, int j, int k, double & T);
  void init();

#if COLLECT_STATS
  olim3d_node_stats & get_node_stats(int i, int j, int k);
  olim3d_node_stats const & get_node_stats(int i, int j, int k) const;
  olim3d_node_stats * _node_stats {nullptr};
#endif
};

template <
  class node, class line_updates, class tri_updates, class tetra_updates,
  class groups>
struct olim3d:
  public abstract_olim3d<
    olim3d<node, line_updates, tri_updates, tetra_updates, groups>,
    node, line_updates, tri_updates, tetra_updates, groups::nneib>,
  public groups
{
  using abstract_olim3d<
    olim3d<node, line_updates, tri_updates, tetra_updates, groups>,
    node, line_updates, tri_updates, tetra_updates,
    groups::nneib>::abstract_olim3d;

  void update_crtp(int i, int j, int k, double & T);
};

template <class groups>
using olim3d_mp0 = olim3d<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates,
  update_rules::mp0_tetra_updates,
  groups>;

template <class groups>
using olim3d_mp1 = olim3d<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates,
  update_rules::mp1_tetra_updates,
  groups>;

template <class groups>
using olim3d_rhr = olim3d<
  node_3d,
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates,
  update_rules::rhr_tetra_updates,
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

template <
  class node, class line_updates, class tri_updates, class tetra_updates>
struct olim3d_hu:
  public abstract_olim3d<
    olim3d_hu<node, line_updates, tri_updates, tetra_updates>,
    node, line_updates, tri_updates, tetra_updates, 26>
{
  using abstract_olim3d<
    olim3d_hu<node, line_updates, tri_updates, tetra_updates>,
    node, line_updates, tri_updates, tetra_updates, 26>::abstract_olim3d;

  void update_crtp(int i, int j, int k, double & T);
};

using olim3d_hu_rhr = olim3d_hu<
  node_3d,
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates,
  update_rules::rhr_tetra_updates>;

using olim3d_hu_mp0 = olim3d_hu<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates,
  update_rules::mp0_tetra_updates>;

using olim3d_hu_mp1 = olim3d_hu<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates,
  update_rules::mp1_tetra_updates>;

#include "olim3d.impl.hpp"

#endif // __OLIM3D_HPP__
