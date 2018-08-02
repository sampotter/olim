#ifndef __HYBRID_IMPL_HPP__
#define __HYBRID_IMPL_HPP__

#include <cassert>

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

  T c = a, fc = f(c), fd, d, dm, df, ds, dd;

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
    ds = df == 0 ? dm : -fb*(a - b)/df;
    dd = sgn(ds) != sgn(dm) || fabs(ds) > fabs(dm) ? dm : ds;
    if (fabs(dd) < tol) {
      dd = tol*sgn(dm)/2;
    }
    d = b + dd;
    fd = f(d);
    if (fd == 0) {
      c = d;
      b = c;
      fc = fd;
      fb = fc;
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
  return {(b + c)/2, hybrid_status::OK};
}

#endif // __HYBRID_IMPL_HPP__
