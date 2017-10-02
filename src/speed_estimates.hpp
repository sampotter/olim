#ifndef __SPEED_ESTIMATES_HPP__
#define __SPEED_ESTIMATES_HPP__

struct rhr_speed_estimate {
  double estimate_speed(double s, double s0) const;
  double estimate_speed(double s, double s0, double s1) const;
  double estimate_speed(double s, double s0, double s1, double s2) const;
};

struct mp0_speed_estimate {
  double estimate_speed(double s, double s0) const;
  double estimate_speed(double s, double s0, double s1) const;
  double estimate_speed(double s, double s0, double s1, double s2) const;
};

#endif // __SPEED_ESTIMATES_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
