#ifndef __HYBRID_IMPL_HPP__
#define __HYBRID_IMPL_HPP__

#include <cassert>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <class F, class T>
std::pair<T, hybrid_status>
hybrid(F const & f, T a, T b, T tol) {
  T c = a, fa = f(a), fb = f(b), fc = f(c), fd;
  T d, dm, df, ds, dd;

  assert(!(fa == 0 || fb == 0 || fc == 0));

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
    a = b = d;
    fa = fb = fd;
    if (sgn(fb) == sgn(fc)) {
      c = a;
      fc = fa;
    }
  }
  return {(b + c)/2, hybrid_status::OK};
}

#endif // __HYBRID_IMPL_HPP__
