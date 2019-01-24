#include "basic_marcher.hpp"

#include <src/config.hpp>

#include <assert.h>
#include <math.h>

#include "common.hpp"

void basic_marcher::update_impl(int lin, double & T) {
  (void) lin;

  double sh = this->s_hat*this->get_h();

  // double T1 = fmin(
  //   nb[0] ? this->nb[0]->get_value() : inf<double>,
  //   nb[2] ? this->nb[2]->get_value() : inf<double>);

  double T1 = fmin(
    nb[0] != -1 ? _U[nb[0]] : inf<double>,
    nb[2] != -1 ? _U[nb[2]] : inf<double>);

  // double T2 = fmin(
  //   nb[1] ? this->nb[1]->get_value() : inf<double>,
  //   nb[3] ? this->nb[3]->get_value() : inf<double>);

  double T2 = fmin(
    nb[1] != -1 ? _U[nb[1]] : inf<double>,
    nb[3] != -1 ? _U[nb[3]] : inf<double>);

  bool T1_inf = isinf(T1), T2_inf = isinf(T2);

  if (!T1_inf && !T2_inf) {
    double diff = T1 - T2, disc = 2*sh*sh - diff*diff;
    T = disc > 0 ? fmin(T, (T1 + T2 + sqrt(disc))/2) : T;
  } else if (T1_inf) {
    T = fmin(T, T2 + sh);
  } else if (T2_inf) {
    T = fmin(T, T1 + sh);
  } else {
    assert(false);
  }
}
