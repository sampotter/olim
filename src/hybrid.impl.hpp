#pragma once

#include <assert.h>

// TODO: almost all of the uses of this in hybrid are to check if the
// sign of two numbers is the same or different: we might be able to
// speed this up by using this function from Stephen Canon on SO:
//
// template <typename T>
// inline bool samesign(T a, T b) {
//   return a*b >= 0;
// }
//
// See here for a few comments:
// https://stackoverflow.com/questions/2922619/how-to-efficiently-compare-the-sign-of-two-floating-point-values-while-handling
//
// This may be worth trying out since it appears that most of the time
// spent *inside* hybrid (i.e., not in calls to f), is spent in `sgn'.
template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

template <class F, class T>
std::pair<T, hybrid_status>
hybrid(F const & f, T a, T b, T tol)
{
  T fa = f(a);
  if (fabs(fa)/fabs(a) <= tol) {
    return {a, hybrid_status::OK};
  }

  T fb = f(b);
  if (fabs(fb)/fabs(b) <= tol) {
    return {b, hybrid_status::OK};
  }

  T c = a, fc = fa, fd, d, dm, df, ds, dd;

  if (sgn(fb) == sgn(fc)) {
    return {0, hybrid_status::DEGENERATE};
  }

  while (true) {
    if (fabs(fc) < fabs(fb)) {
      std::swap(b, c);
      std::swap(fb, fc);
      a = c;
      fa = fc;
    }
    if (fabs(b - c) <= tol) {
      break;
    }
    dm = (c - b)/2;
    df = fa - fb;
    ds = fabs(df) < tol ? dm : -fb*(a - b)/df;
    dd = sgn(ds) != sgn(dm) || fabs(ds) > fabs(dm) ? dm : ds;
    if (fabs(dd) < tol) {
      dd = tol*sgn(dm)/2;
    }
    d = b + dd;
    fd = f(d);
    if (fabs(fd) < tol) {
      b = c = d;
      fb = fc = fd;
      break;
    }
    a = b;
    b = d;
    fa = fb;
    fb = fd;
    if (sgn(fb) == sgn(fc)) {
      c = a;
      fc = fa;
    }
  }

  // NOTE: in G.W. Stewart's Afternotes chapter on this algorithm, he
  // notes that on termination the b and c will satisfy |c - b| <
  // eps. So, we might as well take our final result to be the average
  // of the two in this case.
  return {(b + c)/2, hybrid_status::OK};
}
