#ifndef __OLIM3D_HPP__
#define __OLIM3D_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tetra_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class node, class line_updates, class tri_updates,
          class tetra_updates, class groups>
struct olim3d: public marcher_3d<node>, public line_updates, public tri_updates,
               public tetra_updates
{
  static constexpr bool do_line2_updates =
    groups::group_I || groups::group_II || groups::group_III ||
    groups::group_IV_b || groups::group_V || groups::group_VI_b;

  static constexpr bool do_line3_updates =
    groups::group_V || groups::group_VI_a || groups::group_VI_b;

  static constexpr bool do_tri11_updates =
    groups::group_II || groups::group_III || groups::group_IV_a ||
    groups::group_VI_a;

  static constexpr bool do_tri12_updates =
    groups::group_I || groups::group_II || groups::group_III || groups::group_V;

  static constexpr bool do_tri13_updates = groups::group_V;

  static constexpr bool do_tri22_updates =
    groups::group_I || groups::group_II || groups::group_III ||
    groups::group_IV_b || groups::group_VI_b;

  static constexpr bool do_tri23_updates =
    groups::group_V || groups::group_VI_b;

  using marcher_3d<node>::marcher_3d;
protected:
  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb);
  static int di[26];
  static int dj[26];
  static int dk[26];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
  virtual void update_impl(int i, int j, int k, double & T);
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

struct olim6_groups {
  static constexpr bool group_I = false;
  static constexpr bool group_II = false;
  static constexpr bool group_III = false;
  static constexpr bool group_IV_a = true;
  static constexpr bool group_IV_b = false;
  static constexpr bool group_V = false;
  static constexpr bool group_VI_a = false;
  static constexpr bool group_VI_b = false;
};

using olim6_mp0 = olim3d_mp0<olim6_groups>;
using olim6_mp1 = olim3d_mp1<olim6_groups>;
using olim6_rhr = olim3d_rhr<olim6_groups>;

struct olim18_groups {
  static constexpr bool group_I = true;
  static constexpr bool group_II = true;
  static constexpr bool group_III = true;
  static constexpr bool group_IV_a = true;
  static constexpr bool group_IV_b = true;
  static constexpr bool group_V = false;
  static constexpr bool group_VI_a = false;
  static constexpr bool group_VI_b = false;
};

using olim18_mp0 = olim3d_mp0<olim18_groups>;
using olim18_mp1 = olim3d_mp1<olim18_groups>;
using olim18_rhr = olim3d_rhr<olim18_groups>;

struct olim26_groups {
  static constexpr bool group_I = true;
  static constexpr bool group_II = false;
  static constexpr bool group_III = false;
  static constexpr bool group_IV_a = true;
  static constexpr bool group_IV_b = true;
  static constexpr bool group_V = true;
  static constexpr bool group_VI_a = true;
  static constexpr bool group_VI_b = true;
};

using olim26_mp0 = olim3d_mp0<olim26_groups>;
using olim26_mp1 = olim3d_mp1<olim26_groups>;
using olim26_rhr = olim3d_rhr<olim26_groups>;

#include "olim3d.impl.hpp"

#endif // __OLIM3D_HPP__
