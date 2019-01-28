#pragma once

namespace updates {

template <int d>
struct info;

template <>
struct info<1>
{
  double value {inf<double>};
  vec<double, 1> lambda = {0.5};
  double tol {1e1*eps<double>};
  inline bool inbounds() const {
    return 0 <= lambda[0] && lambda[0] <= 1;
  }
  inline bool finite_lambda() const {
    return isfinite(lambda[0]);
  }
  inline bool in_interior() const {
    return eps<double> <= lambda[0] && lambda[0] <= 1 - eps<double>;
  }
};

template <>
struct info<2>
{
  double value {inf<double>};
  vec2<double> lambda = {1./3, 1./3};
  double tol {1e1*eps<double>};
  inline bool inbounds() const {
    return 0 <= lambda[0] && 0 <= lambda[1] && lambda[0] + lambda[1] <= 1;
  }
  inline bool finite_lambda() const {
    return isfinite(lambda[0]) && isfinite(lambda[1]);
  }
  inline bool in_interior() const {
    return eps<double> <= lambda[0] && eps<double> <= lambda[1] &&
      eps<double> <= 1 - lambda[0] - lambda[1];
  }
};

#ifdef COLLECT_STATS
template <int n>
struct stats
{
  int num_visits {0};
  int count[n];

  stats() { for (int d = 0; d < n; ++d) count[d] = 0; }
};
#endif

}

