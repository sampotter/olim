#ifndef __OLIM8HU_HPP__
#define __OLIM8HU_HPP__

#include "moore_marcher.hpp"
#include "node.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class line_updates,
          class tri_updates>
struct olim8hu: public moore_marcher<node>,
                public line_updates,
                public tri_updates {
  using moore_marcher::moore_marcher;
private:
  virtual void update_impl(int i, int j, double & T);
};

using olim8hu_mp0 = olim8hu<
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates<true>
>;

#include "olim8hu.impl.hpp"

#endif // __OLIM8HU_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
