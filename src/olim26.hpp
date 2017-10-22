#ifndef __OLIM26_HPP__
#define __OLIM26_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tetra_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class node,
          class line_updates,
          class tri_updates,
          class tetra_updates>
struct olim26: public marcher_3d<node>,
               public line_updates,
               public tri_updates,
               public tetra_updates {
  enum dir {
    N, E, U, S, W, D, // degree 1
 // 0  1  2  3  4  5
    UN, UE, US, UW, NE, SE, SW, NW, DN, DE, DS, DW, // degree 2
 // 6   7   8   9   10  11  12  13  14  15  16  17
    UNE, USE, USW, UNW, DNE, DSE, DSW, DNW, // degree 3
 // 18   19   20   21   22   23   24   25
  };

  using marcher_3d<node>::marcher_3d;
protected:
  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb);
  static int di[26];
  static int dj[26];
  static int dk[26];
  static int line1tris[6][8];
  static int line2tris[12][4];
  static int line3tris[8][6];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
  virtual void update_impl(int i, int j, int k, double & T);
};

using olim26_mp0 = olim26<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates<false>,
  update_rules::mp0_tetra_updates
>;

using olim26_rhr = olim26<
  node_3d,
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates<false>,
  update_rules::rhr_tetra_updates
>;

#include "olim26.impl.hpp"

#endif // __OLIM26_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
