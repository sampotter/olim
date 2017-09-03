#ifndef __OLIM_UTIL_HPP__
#define __OLIM_UTIL_HPP__

template <int N>
double rhr(double const * p0, double const * dp, double u0, double u1,
           double s_est, double h);

double rhr_adj(double u0, double u1, double s_est, double h,
               double * lam = nullptr);
double rhr_diag(double u0, double u1, double s_est, double h);
double rhr_3d22(double u0, double u1, double s_est, double h);

#endif // __OLIM_UTIL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
