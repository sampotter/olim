#ifndef __UPDATES_HPP__
#define __UPDATES_HPP__

namespace updates {

template <int d>
struct info;

template <>
struct info<1>
{
  double value {inf<double>};
  double lambda[1] = {0.5};
  double tol {1e1*eps<double>};
  inline bool inbounds() const {
    return 0 <= lambda[0] && lambda[0] <= 1;
  }
  inline bool on_boundary() const {
    return lambda[0] < tol || 1 - lambda[0] < tol;
  }
#ifdef COLLECT_STATS
  bool hierarchical {false};
#endif
};

template <>
struct info<2>
{
  double value {inf<double>};
  double lambda[2] = {1./3, 1./3};
  double tol {1e1*eps<double>};
  inline bool inbounds() const {
    return 0 <= lambda[0] && 0 <= lambda[1] && lambda[0] + lambda[1] <= 1;
  }
  inline bool on_boundary() const {
    return lambda[0] < tol || lambda[1] < tol || lambda[0] + lambda[1] > 1 - tol;
  }
#ifdef COLLECT_STATS
  bool hierarchical {false};
#endif
};

}

#endif // __UPDATES_HPP__
