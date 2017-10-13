#ifndef __OLIM8_MP1_HPP__
#define __OLIM8_MP1_HPP__

#include "olim8.hpp"
#include "olim8hu.hpp"
#include "olim8lut.hpp"

struct olim8_mp1_update_rules {
  double adj1pt(double u0, double s, double s0, double h) const;
  double diag1pt(double u0, double s, double s0, double h) const;
  double adj2pt(double u0, double u1, double s, double s0, double s1, double h) const;
  double diag2pt(double u0, double u1, double s, double s0, double s1, double h) const;
};

using olim8_mp1 = olim8<olim8_mp1_update_rules>;
using olim8hu_mp1 = olim8hu<olim8_mp1_update_rules>;
using olim8lut_mp1 = olim8lut<olim8_mp1_update_rules>;

#endif // __OLIM8_MP1_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
