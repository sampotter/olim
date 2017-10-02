#include "speed_estimates.hpp"

double rhr_speed_estimate::estimate_speed(double s, double s0) const {
  (void) s0;
  return s;
}

double rhr_speed_estimate::estimate_speed(double s, double s0, double s1) const
{
  (void) s0;
  (void) s1;
  return s;
}

double rhr_speed_estimate::estimate_speed(
  double s, double s0, double s1, double s2) const
{
  (void) s0;
  (void) s1;
  (void) s2;
  return s;
}

double mp0_speed_estimate::estimate_speed(double s, double s0) const {
  return (s + s0)/2;
}

double mp0_speed_estimate::estimate_speed(double s, double s0, double s1) const
{
  return (s + (s0 + s1)/2)/2;
}

double mp0_speed_estimate::estimate_speed(
  double s, double s0, double s1, double s2) const
{
  return (s + (s0 + s1 + s2)/3)/2;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
