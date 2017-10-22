#ifndef __OLIM8_HPP__
#define __OLIM8_HPP__

#include "moore_marcher.hpp"
#include "node.hpp"
#include "update_rules.line_updates.hpp"
#include "update_rules.tri_updates.hpp"

template <class line_updates,
          class tri_updates>
struct olim8: public moore_marcher<node>,
              public line_updates,
              public tri_updates {
  using moore_marcher::moore_marcher;
private:
  virtual void update_impl(int i, int j, double & T);
};

using olim8_mp0 = olim8<
  update_rules::mp_line_updates,
  update_rules::mp0_tri_updates<true>
>;

using olim8_mp1 = olim8<
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates<true>
>;

using olim8_rhr = olim8<
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates<true>
>;

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

using olim8hu_mp1 = olim8hu<
  update_rules::mp_line_updates,
  update_rules::mp1_tri_updates<true>
>;

using olim8hu_rhr = olim8hu<
  update_rules::rhr_line_updates,
  update_rules::rhr_tri_updates<true>
>;

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

#include "olim8.impl.hpp"

#endif // __OLIM8_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
