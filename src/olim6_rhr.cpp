#include "olim6_rhr.hpp"

#include "conics.hpp"

void arma_rootfinder::intersect_conics(double const * Q1, double const * Q2,
                                       double * P, int & n) const {
  ::intersect_conics(Q1, Q2, P, n);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
