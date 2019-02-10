#include "basic_marcher.hpp"

#include <src/config.hpp>

#include <assert.h>
#include <math.h>

#include "common.hpp"

void basic_marcher::update_impl(int lin, int const * nb, int parent, double & U) {
  (void) lin;
  (void) parent;

  double sh = this->_s[lin]*this->get_h();

  double U1 = fmin(
    nb[0] != -1 ? _U[nb[0]] : inf<double>,
    nb[2] != -1 ? _U[nb[2]] : inf<double>);

  double U2 = fmin(
    nb[1] != -1 ? _U[nb[1]] : inf<double>,
    nb[3] != -1 ? _U[nb[3]] : inf<double>);

  bool U1_inf = isinf(U1), U2_inf = isinf(U2);

  if (!U1_inf && !U2_inf) {
    double diff = U1 - U2, disc = 2*sh*sh - diff*diff;
    U = disc > 0 ? fmin(U, (U1 + U2 + sqrt(disc))/2) : U;
  } else if (U1_inf) {
    U = fmin(U, U2 + sh);
  } else if (U2_inf) {
    U = fmin(U, U1 + sh);
  } else {
    assert(false);
  }
}
