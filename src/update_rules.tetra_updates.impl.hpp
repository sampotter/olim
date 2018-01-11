#ifndef __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__
#define __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.defs.hpp"
#include "common.macros.hpp"
#include "numopt.hpp"
#include "olim_util.hpp"
#include "update_rules.tetra_updates.util.hpp"

template <class speed_est, char degree>
template <char p0, char p1, char p2>
double update_rules::tetra_updates<speed_est, degree>::tetra(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  ffvec<p0>, ffvec<p1>, ffvec<p2>) const
{
  return tetra_impl<p0, p1, p2>(
    u0, u1, u2, s, s0, s1, s2, h,
    ffvec<p0> {}, ffvec<p1> {}, ffvec<p2> {},
    std::integral_constant<char, degree> {});
}

template <class speed_est, char degree>
double update_rules::tetra_updates<speed_est, degree>::tetra(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  return tetra_impl(
    p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h,
    std::integral_constant<char, degree> {});
}

template <class speed_est, char degree>
template <char p0, char p1, char p2>
double update_rules::tetra_updates<speed_est, degree>::tetra_impl(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  ffvec<p0>, ffvec<p1>, ffvec<p2>,
  std::integral_constant<char, 0>) const
{
  double p0_vec[3] = {component(p0, 0), component(p0, 1), component(p0, 2)};
  double p1_vec[3] = {component(p1, 0), component(p1, 1), component(p1, 2)};
  double p2_vec[3] = {component(p2, 0), component(p2, 1), component(p2, 2)};
  return tetra(p0_vec, p1_vec, p2_vec, u0, u1, u2, s, s0, s1, s2, h);
}

template <class speed_est, char degree>
template <char p0, char p1, char p2>
double update_rules::tetra_updates<speed_est, degree>::tetra_impl(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  ffvec<p0>, ffvec<p1>, ffvec<p2>,
  std::integral_constant<char, 1>) const
{
  double p0_vec[3] = {component(p0, 0), component(p0, 1), component(p0, 2)};
  double p1_vec[3] = {component(p1, 0), component(p1, 1), component(p1, 2)};
  double p2_vec[3] = {component(p2, 0), component(p2, 1), component(p2, 2)};
  return tetra(p0_vec, p1_vec, p2_vec, u0, u1, u2, s, s0, s1, s2, h);
}

template <class speed_est, char degree>
double update_rules::tetra_updates<speed_est, degree>::tetra_impl(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  std::integral_constant<char, 0>) const
{
#ifdef USE_ARMADILLO
  using namespace arma;
  using namespace numopt;

  vec const P0(const_cast<double *>(p0), 3, false, true);
  vec const P1(const_cast<double *>(p1), 3, false, true);
  vec const P2(const_cast<double *>(p2), 3, false, true);

  double sh = this->s_hat(s, s0, s1, s2)*h;
  mat const dP = arma::join_horiz(P1 - P0, P2 - P0);
  vec const du = {u1 - u0, u2 - u0};

  field_t const F0 = [&] (vec const & x) -> double {
    return (1 - sum(x))*u0 + x(0)*u1 + x(1)*u2 +
      sh*norm((1 - sum(x))*P0 + x(0)*P1 + x(1)*P2);
  };

  grad_t const dF0 = [&] (vec const & x) -> vec {
    vec const P = (1 - sum(x))*P0 + x(0)*P1 + x(1)*P2;
    return du + sh*dP.t()*P/norm(P);
  };

  hess_t const d2F0 = [&] (vec const & x) -> mat {
    vec const P = (1 - sum(x))*P0 + x(0)*P1 + x(1)*P2;
    return sh*dP.t()*(eye(3, 3) - P*P.t()/dot(P, P))*dP/norm(P);
  };

  mat const A = {{1, 0}, {0, 1}, {-1, -1}};
  vec const b = {0, 0, -1};
  vec const x0 = {1./3, 1./3};

  bool error;
  vec const xopt = sqp(F0, dF0, d2F0, A, b, x0, &error);
  assert(!error);

  return F0(xopt);
#else
  double u[3] = {u0, u1, u2};
  double s_hat = s;
  double s_[3] = {s0, s1, s2};
  double p[3][3];
  memcpy((void *) p[0], (void *) p0, 3*sizeof(double));
  memcpy((void *) p[1], (void *) p1, 3*sizeof(double));
  memcpy((void *) p[2], (void *) p2, 3*sizeof(double));
  F0<2> func {h, this->theta()};
  func.set_args(u, s_hat, s_, p);
  double lambda[2], F0;
  bool error;
  numopt::sqp_baryplex(&func, lambda, &error);
  assert(!error);
  func.set_lambda(lambda); // TODO: maybe unnecessary
  func.eval(F0);
  return F0;
#endif
}

template <class speed_est, char degree>
double update_rules::tetra_updates<speed_est, degree>::tetra_impl(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  std::integral_constant<char, 1>) const
{
#ifdef USE_ARMADILLO
  using namespace arma;
  using namespace numopt;

  vec const P0(const_cast<double *>(p0), 3, false, true);
  vec const P1(const_cast<double *>(p1), 3, false, true);
  vec const P2(const_cast<double *>(p2), 3, false, true);

  mat const dP = arma::join_horiz(P1 - P0, P2 - P0);
  vec const du = {u1 - u0, u2 - u0};
  vec const ds = {s1 - s0, s2 - s0};

  auto stheta = [&] (vec const & x) -> double {
    return (1 - this->theta())*s +
      this->theta()*((1 - sum(x))*s0 + x(0)*s1 + x(1)*s2);
  };

  field_t const F1 = [&] (vec const & x) -> double {
    double const lam0 = 1 - sum(x);
    vec const P = lam0*P0 + x(0)*P1 + x(1)*P2;
    double const u = lam0*u0 + x(0)*u1 + x(1)*u2;
    return u + h*stheta(x)*norm(P);
  };

  grad_t const dF1 = [&] (vec const & x) -> vec {
    vec const P = (1 - sum(x))*P0 + x(0)*P1 + x(1)*P2;
    double const l = norm(P);
    return du + h*(this->theta()*l*l*ds + stheta(x)*dP.t()*P)/l;
  };

  hess_t const d2F1 = [&] (vec const & x) -> mat {
    vec const P = (1 - sum(x))*P0 + x(0)*P1 + x(1)*P2;
    double const l = norm(P);
    return h*(this->theta()*(dP.t()*P*ds.t() + ds*P.t()*dP) +
              stheta(x)*dP.t()*(eye(3, 3) - P*P.t()/dot(P, P))*dP)/l;
  };

  mat const A = {{1, 0}, {0, 1}, {-1, -1}};
  vec const b = {0, 0, -1};
  vec const x0 = {1./3, 1./3};

  bool error;
  vec const xopt = sqp(F1, dF1, d2F1, A, b, x0, &error);
  assert(!error);

  return F1(xopt);
#else
  double u[3] = {u0, u1, u2};
  double s_hat = s;
  double s_[3] = {s0, s1, s2};
  double p[3][3];
  memcpy((void *) p[0], (void *) p0, 3*sizeof(double));
  memcpy((void *) p[1], (void *) p1, 3*sizeof(double));
  memcpy((void *) p[2], (void *) p2, 3*sizeof(double));
  F1<2> func {h, this->theta()};
  func.set_args(u, s_hat, s_, p);
  double lambda[2], F1;
  bool error;
  numopt::sqp_baryplex(&func, lambda, &error);
  assert(!error);
  func.set_lambda(lambda); // TODO: maybe unnecessary
  func.eval(F1);
  return F1;
#endif
}

#endif // __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__
