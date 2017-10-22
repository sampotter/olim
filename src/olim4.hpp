#ifndef __OLIM4_HPP__
#define __OLIM4_HPP__

#include "neumann_marcher.hpp"
#include "node.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class line_updates,
          class tri_updates>
struct olim4: public neumann_marcher<node>,
              public line_updates,
              public tri_updates {
  using neumann_marcher::neumann_marcher;
private:
  virtual void update_impl(int i, int j, double & T);
};

using olim4_mp0 = olim4<
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates<true>
>;

using olim4_mp1 = olim4<
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates<true>
>;

using olim4_rhr = olim4<
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates<true>
>;

#include "olim4.impl.hpp"

#endif // __OLIM4_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
