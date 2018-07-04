#ifndef __OLIM_HPP__
#define __OLIM_HPP__

#include <type_traits>

#include "marcher.hpp"
#include "node.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class node, class line_updates, class tri_updates, bool adj_updates,
          bool diag_updates>
struct olim: public marcher<
               olim<
                 node, line_updates, tri_updates, adj_updates,
                 diag_updates>,
               node>,
             public line_updates,
             public tri_updates {
  using marcher<
    olim<
      node, line_updates, tri_updates, adj_updates,
      diag_updates>,
    node>::marcher;
  static_assert(adj_updates || diag_updates, "error");
  static constexpr int nneib = diag_updates ? 8 : 4;
EIKONAL_PRIVATE:
  virtual void update_impl(node * n, node ** nb, double & T);
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
