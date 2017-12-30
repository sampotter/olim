#ifndef __OLIM_HPP__
#define __OLIM_HPP__

#include "neumann_marcher.hpp"
#include "node.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class node, class line_updates, class tri_updates, bool adj_updates,
          bool diag_updates>
struct olim: public neumann_marcher<node>, public line_updates,
             public tri_updates {
  static_assert(adj_updates || diag_updates, "error");
  using neumann_marcher<node>::neumann_marcher;
protected:
  static int di[8];
  static int dj[8];
private:
  virtual void update_impl(int i, int j, double & T);
};

using olim4_mp0 = olim<
  node,
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates,
  true,
  false>;

using olim4_mp1 = olim<
  node,
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates,
  true,
  false>;

using olim4_rhr = olim<
  node,
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates,
  true,
  false>;

using olim8_mp0 = olim<
  node,
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates,
  true,
  true>;

using olim8_mp1 = olim<
  node,
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates,
  true,
  true>;

using olim8_rhr = olim<
  node,
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates,
  true,
  true>;

#include "olim.impl.hpp"

#endif // __OLIM_HPP__
