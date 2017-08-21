#include "rpoly.hpp"

extern "C"
void
rpoly_(double op[], int * degree, double zeror[], double zeroi[], int * fail);

bool rpoly(double * const coefs, int degree, std::complex<double> * roots) {
  double * tmp = new double[3*degree + 1];
  for (int i = 0; i <= degree; ++i) {
    tmp[i] = coefs[degree - i];
  }
  int fail = 0;
  double * zeror = &tmp[degree + 1];
  double * zeroi = &zeror[degree];
  rpoly_(tmp, &degree, &tmp[degree + 1], &tmp[2*degree + 1], &fail);
  for (int i = 0; i < degree; ++i) {
    roots[i] = {zeror[i], zeroi[i]};
  }
  return static_cast<bool>(fail);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
