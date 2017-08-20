#ifndef __QROOTS_HPP__
#define __QROOTS_HPP__

void qroots(double * a, double * roots, double l = 0, double r = 1);

// TODO: for testing only
int sigma(double ** polys, double x);
int oldsigma(double ** polys, double x);
int sturm(double ** polys, double l, double r);

#endif // __QROOTS_HPP__
