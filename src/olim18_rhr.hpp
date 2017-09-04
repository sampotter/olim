#ifndef __OLIM18_RHR_HPP__
#define __OLIM18_RHR_HPP__

#include "conics.hpp"
#include "olim18.hpp"
#include "node_3d.hpp"

template <class rootfinder>
struct olim18_rhr_update_rules: rootfinder {
  double line1(double u0, double s, double s0, double h) const;
  double line2(double u0, double s, double s0, double h) const;
  double tri12(double u0, double u1, double s, double s0, double s1,
               double h) const;
  double tri22(double u0, double u1, double s, double s0, double s1,
               double h) const;
  double tetra122(double u0, double u1, double u2, double s, double s0,
                  double s1, double s2, double h) const;
  double tetra222(double u0, double u1, double u2, double s, double s0,
                  double s1, double s2, double h) const;
};

using olim18_rhr_arma = olim18<node_3d, olim18_rhr_update_rules<arma_rootfinder>>;

#include "olim18_rhr.impl.hpp"

#endif // __OLIM18_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
