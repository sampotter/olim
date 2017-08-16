#ifndef __OLIM8_MP1_HPP__
#define __OLIM8_MP1_HPP__

#include "olim8.hpp"

template <class Rootfinder>
struct olim8_mp1_update_rules: public Rootfinder {
  double adj1pt(double u0, double s, double s0, double h) const;
  double diag1pt(double u0, double s, double s0, double h) const;
  double adj2pt(double u0, double u1, double s, double s0, double s1, double h) const;
  double diag2pt(double u0, double u1, double s, double s0, double s1, double h) const;
};

struct bsearch_rootfinder {
  void find_roots(double * a, double * roots) const;
};

struct gsl_rootfinder {
  void find_roots(double * a, double * roots) const;
};

using olim8_mp1_bsearch = olim8<olim8_mp1_update_rules<bsearch_rootfinder>>;
using olim8_mp1_gsl = olim8<olim8_mp1_update_rules<gsl_rootfinder>>;

#include "olim8_mp1.impl.hpp"

#endif // __OLIM8_MP1_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
