#ifndef __QROOTS_HPP__
#define __QROOTS_HPP__

void qroots(double const * const a, double * roots, double l = 0, double r = 1);

// TODO: for testing only
int sigma(double const * const * polys, double x);
int oldsigma(double const * const * polys, double x);
int sturm(double const * const * polys, double l, double r);

#endif // __QROOTS_HPP__
