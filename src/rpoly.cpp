#include "rpoly.hpp"

#include <cstdio>

extern "C"
void
rpoly_(double op[], int * degree, double zeror[], double zeroi[], int * fail);

bool rpoly(double const * coefs, int degree, std::complex<double> * roots) {
  // TODO: optimize using a static variable and realloc
  double * tmp = new double[3*degree + 1];
  for (int i = 0; i <= degree; ++i) {
    tmp[i] = coefs[degree - i];
  }
  int fail = 0;
  double * zeror = &tmp[degree + 1];
  double * zeroi = &tmp[2*degree + 1];
  rpoly_(tmp, &degree, zeror, zeroi, &fail);
  for (int i = 0; i < degree; ++i) {
    roots[i] = {zeror[i], zeroi[i]};
  }
  delete[] tmp;
  return static_cast<bool>(fail);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
