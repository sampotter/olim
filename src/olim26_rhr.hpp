#ifndef __OLIM26_RHR_HPP__
#define __OLIM26_RHR_HPP__

#include "conics.hpp"
#include "olim26.hpp"
#include "node_3d.hpp"

template <class rootfinder>
struct olim26_rhr_update_rules: rootfinder {
  double line1(double u0, double s, double s0, double h) const;
  double line2(double u0, double s, double s0, double h) const;
  double line3(double u0, double s, double s0, double h) const;
  double tri12(double u0, double u1, double s, double s0, double s1,
               double h) const;
  double tri13(double u0, double u1, double s, double s0, double s1,
               double h) const;
  double tri23(double u0, double u1, double s, double s0, double s1,
               double h) const;
  double tetra123(double u0, double u1, double u2, double s, double s0,
                  double s1, double s2, double h) const;
};

using olim26_rhr_arma = olim26<node_3d, olim26_rhr_update_rules<arma_rootfinder>>;

#include "olim26_rhr.impl.hpp"

#endif // __OLIM26_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
