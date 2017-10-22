#ifndef __OLIM8LUT_HPP__
#define __OLIM8LUT_HPP__

#include "moore_marcher.hpp"
#include "node.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class line_updates,
          class tri_updates>
struct olim8lut: public moore_marcher<node>,
                 public line_updates,
                 public tri_updates {
  using moore_marcher::moore_marcher;
private:
  virtual void update_impl(int i, int j, double & T);
};

using olim8lut_mp0 = olim8lut<
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates<true>
>;

using olim8lut_mp1 = olim8lut<
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates<true>
>;

using olim8lut_rhr = olim8lut<
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates<true>
>;

#include "olim8lut.impl.hpp"

#endif // __OLIM8LUT_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
