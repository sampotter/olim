#include "test.hpp"

#include <cstdio>

#include "conics.hpp"

void intersect_conics_works_1() {
  double P[12], Q[12] = {
    11.98801141019249,
    6.878959176895175,
    33.92041333725862,
    10.62559748131321,
    15.73464971461059,
    33.92041333725862,
    16.29759019485503,
    4.533595214828539,
    33.92041333725862,
    6.316018696650761,
    18.08001367667718,
    33.92041333725862,
  };
  double Q1[6] = {
    3.955,
    3.955,
    0.955,
    -3.955,
    -1.955,
    0.9775
  };
  double Q2[6] = {
    0.68,
    3.68,
    3.68,
    -1.68,
    -3.68,
    0.84
  };
  int n;
  intersect_conics(Q1, Q2, P, n);
  for (int i = 0, j = 1, k = 2; i < 12; i += 3, j += 3, k += 3) {
    IS_APPROX_EQUAL(P[i]/P[k], Q[i]/Q[k]);
    IS_APPROX_EQUAL(P[j]/P[k], Q[j]/Q[k]);
  }
}

int main() {
  intersect_conics_works_1();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
