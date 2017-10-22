#ifndef __OLIM8_RHR_HPP__
#define __OLIM8_RHR_HPP__

#include "olim8.hpp"
#include "olim8hu.hpp"
#include "olim8lut.hpp"

#include "update_rules.line_updates.hpp"

struct olim8_rhr_update_rules: public update_rules::rhr_line_updates {
  double tri11(double u0, double u1, double s, double s0, double s1, double h)
    const;
  double tri12(double u0, double u1, double s, double s0, double s1, double h)
    const;
};

using olim8_rhr = olim8<olim8_rhr_update_rules>;
using olim8hu_rhr = olim8hu<olim8_rhr_update_rules>;
using olim8lut_rhr = olim8lut<olim8_rhr_update_rules>;

#endif // __OLIM8_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
