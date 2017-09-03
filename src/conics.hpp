#ifndef __CONICS_HPP__
#define __CONICS_HPP__

void intersect_conics(double const * Q1, double const * Q2, double * P, int & n);

struct arma_rootfinder {
  void intersect_conics(double const * Q1, double const * Q2, double * P,
                        int & n) const;
};

#endif // __CONICS_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
