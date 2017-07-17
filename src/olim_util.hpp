#ifndef __OLIM_UTIL_HPP__
#define __OLIM_UTIL_HPP__

double rhr_adj(double u0, double u1, double s_est, double h,
               double * lam = nullptr);
double rhr_diag(double u0, double u1, double s_est, double h);

double mp1_adj(double u0, double u1, double sbar0, double sbar1, double h);
double mp1_diag(double u0, double u1, double sbar0, double sbar1, double h);

#endif // __OLIM_UTIL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
