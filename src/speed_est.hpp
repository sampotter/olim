#ifndef __SPEED_ESTS_HPP__
#define __SPEED_ESTS_HPP__

struct speed_est {
  speed_est(double theta): theta_ {theta} {}

  inline double s_hat(double s, double s0) const {
    return s + theta_*(s0 - s);
  }

  inline double s_hat(double s, double s0, double s1) const {
    return s + theta_*((s0 + s1)/2 - s);
  }

  inline double s_hat(double s, double s0, double s1, double s2) const {
    return s + theta_*((s0 + s1 + s2)/3 - s);
  }

  inline double theta() const { return theta_; }

private:
  double theta_;
};

/**
 * This specialization exists to simplify the arithmetic used to
 * compute s_hat when theta == 0.0.
 */
struct rhr_speed_est: public speed_est {
  rhr_speed_est(): speed_est {0.0} {}

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

struct mp_speed_est: public speed_est {
  mp_speed_est(): speed_est {0.5} {}
};

#endif // __SPEED_ESTS_HPP__
