#ifndef __OLIM8_RHR_HPP__
#define __OLIM8_RHR_HPP__

#include "olim8.hpp"
#include "olim8fast.hpp"

struct olim8_rhr_update_rules {
  double adj1pt(double u0, double s, double h) const;
  double adj2pt(double u0, double u1, double s, double h) const;
  double diag1pt(double u0, double s, double h) const;
  double diag2pt(double u0, double u1, double s, double h) const;
};

template <class olim = olim8>
using olim8_rhr = olim<olim8_rhr_update_rules>;

#endif // __OLIM8_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
