#ifndef __CONICS_HPP__
#define __CONICS_HPP__

#ifdef EIKONAL_DEBUG
#    include <armadillo>
#endif

bool intersect_conics(double const * Q1, double const * Q2, double * P, int & n);

#ifdef EIKONAL_DEBUG
arma::mat sym_adjoint(arma::mat const & A);
bool split_deg_conic(arma::mat const & A, arma::vec & m, arma::vec & l);
#endif

struct arma_rootfinder {
  bool intersect_conics(double const * Q1, double const * Q2, double * P,
                        int & n) const;
};

#endif // __CONICS_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
