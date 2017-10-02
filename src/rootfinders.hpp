#ifndef __ROOTFINDERS_HPP__
#define __ROOTFINDERS_HPP__

struct bsearch_rootfinder {
  void find_roots(double const * a, double * roots) const;
};

struct gsl_rootfinder {
  void find_roots(double const * a, double * roots) const;
};

#endif // __ROOTFINDERS_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
