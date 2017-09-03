#ifndef __OLIM6_RHR_HPP__
#define __OLIM6_RHR_HPP__

#include "conics.hpp"
#include "olim6.hpp"

template <class rootfinder>
struct olim6_rhr_update_rules: rootfinder {
  double adj1pt(double u0, double s, double s0, double h) const;
  double adj2pt(double u0, double u1, double s, double s0, double s1,
                double h) const;
  double adj3pt(double u0, double u1, double u2, double s, double s0, double s1,
                double s2, double h) const;
};

using olim6_rhr_arma = olim6<olim6_rhr_update_rules<arma_rootfinder>>;

#include "olim6_rhr.impl.hpp"

#endif // __OLIM6_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
