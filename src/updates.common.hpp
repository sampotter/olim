#ifndef __UPDATES_HPP__
#define __UPDATES_HPP__

namespace updates {

template <int d>
struct info {
  double value {INF(double)};
  double lambda[d];
#ifdef COLLECT_STATS
  inline bool is_degenerate() const {
    bool deg = false;
    double lam0 = 1;
    for (int i = 0; i < d; ++i) {
      deg = deg || lambda[i] < EPS(double);
      lam0 -= lambda[i];
    }
    return deg || lam0 > 1 - EPS(double);
  }
  bool hierarchical {false};
#endif
};

}

#endif // __UPDATES_HPP__
