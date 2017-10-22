#ifndef __OLIM8_MP0_HPP__
#define __OLIM8_MP0_HPP__

#include "olim8.hpp"
#include "olim8hu.hpp"
#include "olim8lut.hpp"

#include "update_rules.line_updates.hpp"

struct olim8_mp0_update_rules {
  double tri11(double u0, double u1, double s, double s0, double s1, double h) const;
  double tri12(double u0, double u1, double s, double s0, double s1, double h) const;
};

using olim8_mp0 = olim8<update_rules::mp_line_updates, olim8_mp0_update_rules>;
using olim8hu_mp0 = olim8hu<update_rules::mp_line_updates, olim8_mp0_update_rules>;
using olim8lut_mp0 = olim8lut<update_rules::mp_line_updates, olim8_mp0_update_rules>;

#endif // __OLIM8_MP0_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
