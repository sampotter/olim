#ifndef __OLIM18_HPP__
#define __OLIM18_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tetra_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class node,
          class line_updates,
          class tri_updates,
          class tetra_updates>
struct olim18_rect: public marcher_3d<node>,
                    public line_updates,
                    public tri_updates,
                    public tetra_updates {
  using marcher_3d<node>::marcher_3d;
protected:
  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb);
  static int di[18];
  static int dj[18];
  static int dk[18];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
  virtual void update_impl(int i, int j, int k, double & T);
  void do_tri12_updates(
    abstract_node const * const * nb, int const * dirs,
    double const * s_, double s, double h, double & T) const;
  void do_tri22_updates(
    abstract_node const * const * nb, int const * dirs,
    double const * s_, double s, double h, double & T) const;
};

using olim18_mp0 = olim18_rect<
  node_3d,
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates<false>,
  update_rules::mp0_tetra_updates
>;

using olim18_rhr = olim18_rect<
  node_3d,
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates<false>,
  update_rules::rhr_tetra_updates
>;

#include "olim18.impl.hpp"

#endif // __OLIM18_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
