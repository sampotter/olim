#ifndef __OLIM_UTIL_HPP__
#define __OLIM_UTIL_HPP__

#ifdef EIKONAL_DEBUG
void check_params(double u0, double u1, double s, double h);
void check_params(double u0, double u1, double s, double s0, double s1,
                  double h);
void check_params(double u0, double u1, double u2, double s, double h);
void check_params(double u0, double u1, double u2, double s,
                  double s0, double s1, double s2, double h);
#endif

template <int N>
double rhr(double const * p0, double const * dp, double u0, double u1,
           double s_est, double h);

double rhr_adj(double u0, double u1, double s_est, double h,
               double * lam = nullptr);
double rhr_diag(double u0, double u1, double s_est, double h);

#endif // __OLIM_UTIL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
