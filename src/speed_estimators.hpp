#ifndef __SPEED_ESTIMATORS_HPP__
#define __SPEED_ESTIMATORS_HPP__

struct rhr_speed_estimator {
  inline double s_hat(double s, double s0) const {
    (void) s0;
    return s;
  }

  inline double s_hat(double s, double s0, double s1) const {
    (void) s0;
    (void) s1;
    return s;
  }

  inline double s_hat(double s, double s0, double s1, double s2) const {
    (void) s0;
    (void) s1;
    (void) s2;
    return s;
  }
};

struct mp_speed_estimator {
  inline double s_hat(double s, double s0) const {
    return (s + s0)/2;
  }

  inline double s_hat(double s, double s0, double s1) const {
    return (s + (s0 + s1)/2)/2;
  }

  inline double s_hat(double s, double s0, double s1, double s2) const {
    return (s + (s0 + s1 + s2)/3)/2;
  }
};

#endif // __SPEED_ESTIMATORS_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
