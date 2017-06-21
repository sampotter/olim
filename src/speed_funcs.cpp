#include "speed_funcs.hpp"

#include <cmath>

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

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
