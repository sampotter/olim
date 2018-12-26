#ifndef __UPDATES_HPP__
#define __UPDATES_HPP__

namespace updates {

template <int d>
struct info {
  double value {INF(double)};
  double lambda[d];
  inline bool is_degenerate() const {
    bool deg = false;
    double lam0 = 1;
    for (int i = 0; i < d; ++i) {
      deg = deg || lambda[i] < EPS(double);
      lam0 -= lambda[i];
    }
    return std::isinf(value) || deg || lam0 > 1 - EPS(double);
  }
#ifdef COLLECT_STATS
  bool hierarchical {false};
#endif
};

template <>
struct info<2> {
  double value {INF(double)};
  double lambda[2] = {1./3, 1./3};
  inline bool is_degenerate() const {
    return std::isinf(value) ||
      lambda[0] < 0 || lambda[1] < 0 || lambda[0] + lambda[1] > 1;
  }
#ifdef COLLECT_STATS
  bool hierarchical {false};
#endif
};

}

#endif // __UPDATES_HPP__
