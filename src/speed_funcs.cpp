#include "speed_funcs.hpp"

#include <cmath>

double default_speed_func(double x, double y) {
  (void) x;
  (void) y;
  return 1.0;
}

double default_speed_func_soln(double x, double y) {
  return std::sqrt(x*x + y*y);
}

double sf1(double x, double y) {
  double a = (x + y)/2;
  double b = std::cos(a);
  double c = std::sin(b)/2;
  double vx = (x + b)*(1 - c);
  double vy = y - b*c;
  return std::sqrt(vx*vx + vy*vy);
}

double sf1_soln(double x, double y) {
  double a = x + std::cos((x + y)/2);
  return (a*a + y*y)/2;
}

double sf2(double x, double y) {
  double x2 = x*x, y2 = y*y;
  double denom = std::pow(x2 + y2, 4.5);
  double sum = x + y;
  double diff = x - y;
  double x4 = x2*x2, y4 = y2*y2;
  double vx = -4*x*y2*sum*diff*(x4 - 11*x2*y2 + 2*y4)/denom;
  double vy = 4*x2*y*sum*diff*(2*x4 - 11*x2*y2 + y4)/denom;
  return std::sqrt(vx*vx + vy*vy);
}

double sf2_soln(double x, double y) {
  double x2 = x*x, y2 = y*y;
  return 4*x2*y2*std::pow(x2 + y2, 2)/std::pow(x2 - y2, 3.5) + 1;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
