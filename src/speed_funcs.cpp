#include "speed_funcs.hpp"

#include <cmath>

#include "common.hpp"

using std::cos;
using std::fabs;
using std::pow;
using std::sin;
using std::sqrt;

static double r(double x, double y) {
  return sqrt(x*x + y*y);
}

static double r(double x, double y, double z) {
  return sqrt(x*x + y*y + z*z);
}

double default_speed_func(double x, double y) {
  (void) x;
  (void) y;
  return 1.0;
}

double default_speed_func(double x, double y, double z) {
  (void) x;
  (void) y;
  (void) z;
  return 1.0;
}

double default_speed_func_soln(double x, double y) {
  return r(x, y);
}

double default_speed_func_soln(double x, double y, double z) {
  return r(x, y, z);
}

double s1(double x, double y) {
  return 1 - sin(r(x, y));
}

double f1(double x, double y) {
  return cos(r(x, y)) + r(x, y) - 1;
}

double s1(double x, double y, double z) {
  return 1 - sin(r(x, y, z));
}

double f1(double x, double y, double z) {
  return cos(r(x, y, z)) + r(x, y, z) - 1;
}

double s3(double x, double y) {
  return sqrt(x*x + y*y);
}

double f3(double x, double y) {
  return (x*x + y*y)/2;
}

double s5(double x, double y) {
  return sqrt(pow(x, 18) + pow(y, 18));
}

double f5(double x, double y) {
  return (pow(x, 10) + pow(y, 10))/10;
}

double s7(double x, double y) {
  x = x - 0.900367222589747;
  double aux0 = (x + y)/2, aux1 = x + cos(aux0), aux2 = sin(aux0)/2;
  double dx = aux1*(1 - aux2), dy = y - aux1*aux2;
  return sqrt(dx*dx + dy*dy);
}

double f7(double x, double y) {
  x = x - 0.900367222589747;
  return (y*y + pow(x + cos((x + y)/2), 2))/2;
}
