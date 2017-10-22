#ifndef __OLIM6_HPP_HPP__
#define __OLIM6_HPP_HPP__

#include "neumann_marcher_3d.hpp"
#include "node_3d.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tetra_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class line_updates,
          class tri_updates,
          class tetra_updates>
struct olim6: public neumann_marcher_3d<node_3d>,
              public line_updates,
              public tri_updates,
              public tetra_updates {
  enum neighbor {N, E, U, S, W, D};
  
  using neumann_marcher_3d::neumann_marcher_3d;
private:
  virtual void update_impl(int i, int j, int k, double & T);
};

using olim6_mp0 = olim6<
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates,
  update_rules::mp0_tetra_updates
>;

using olim6_rhr = olim6<
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates,
  update_rules::rhr_tetra_updates
>;

#include "olim6.impl.hpp"

#endif // __OLIM6_HPP_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
