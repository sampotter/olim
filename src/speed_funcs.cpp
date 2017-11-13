#include "speed_funcs.hpp"

#include <cmath>

#include "common.defs.hpp"

static double r(double x, double y) {
  return std::sqrt(x*x + y*y);
}

static double r(double x, double y, double z) {
  return std::sqrt(x*x + y*y + z*z);
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

double default_speed_func_soln_3d(double x, double y, double z) {
  return r(x, y, z);
}

double s1(double x, double y) {
  return 1 - std::sin(r(x, y));
}

double f1(double x, double y) {
  return std::cos(r(x, y)) + r(x, y) - 1;
}

double s1(double x, double y, double z) {
  return 1 - std::sin(r(x, y, z));
}

double f1(double x, double y, double z) {
  return std::cos(r(x, y, z)) + r(x, y, z) - 1;
}

double s2(double x, double y) {
  return std::fabs(x + y);
}

double f2(double x, double y) {
  return (x + y)*(x + y)/(2*sqrt2);
}

double s3(double x, double y) {
  return std::sqrt(x*x + y*y);
}

double f3(double x, double y) {
  return (x*x + y*y)/2;
}

double s4(double x, double y) {
  double x_sq = x*x;
  double y_sq = y*y;
  return 2*std::sqrt(
    std::pow(x_sq - y_sq, 2)*(x_sq*x_sq + 14*x_sq*y_sq + y_sq*y_sq)/
    std::pow(x_sq + y_sq, 3));
}

double f4(double x, double y) {
  double x_sq = x*x;
  double y_sq = y*y;
  return -3*x_sq + y_sq + 4*x_sq*x_sq/(x_sq + y_sq);
}

double s5(double x, double y) {
  return std::sqrt(std::pow(x, 18) + std::pow(y, 18));
}

double f5(double x, double y) {
  return (std::pow(x, 10) + std::pow(y, 10))/10;
}

double s6(double x, double y) {
  return 4*std::sqrt(std::pow(x - y, 2)*std::pow(x + y, 2)*(x*x + y*y));
}

double f6(double x, double y) {
  return std::pow(x*x - y*y, 2);
}

double s7(double x, double y) {
  x = x - 0.900367222589747;
  double aux0 = (x + y)/2, aux1 = x + std::cos(aux0), aux2 = std::sin(aux0)/2;
  double dx = aux1*(1 - aux2), dy = y - aux1*aux2;
  return std::sqrt(dx*dx + dy*dy);
}

double f7(double x, double y) {
  x = x - 0.900367222589747;
  return (y*y + std::pow(x + cos((x + y)/2), 2))/2;
}
