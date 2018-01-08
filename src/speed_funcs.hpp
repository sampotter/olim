#ifndef __SPEED_FUNCS_HPP__
#define __SPEED_FUNCS_HPP__

#include <array>

#include "typedefs.h"

struct no_speed_func_t {};
constexpr no_speed_func_t no_speed_func {};

double default_speed_func(double x, double y);
double default_speed_func_soln(double x, double y);

double default_speed_func(double x, double y, double z);
double default_speed_func_soln(double x, double y, double z);

double s1(double x, double y);
double f1(double x, double y);

double s1(double x, double y, double z);
double f1(double x, double y, double z);

double s2(double x, double y);
double f2(double x, double y);

double s3(double x, double y);
double f3(double x, double y);

double s4(double x, double y);
double f4(double x, double y);

double s5(double x, double y);
double f5(double x, double y);

double s6(double x, double y);
double f6(double x, double y);

// Masha's speed function
double s7(double x, double y);
double f7(double x, double y);

static std::array<speed_func, 8> speed_funcs {{
  default_speed_func, s1, s2, s3, s4, s5, s6, s7
}};

static std::array<double(*)(double, double), 8> speed_func_solns {{
  default_speed_func_soln, f1, f2, f3, f4, f5, f6, f7
}};

static std::array<double(*)(double, double, double), 2> speed_funcs_3d {{
  default_speed_func, s1
}};

static std::array<double(*)(double, double, double), 2> speed_func_solns_3d {{
  default_speed_func_soln, f1
}};

#endif // __SPEED_FUNCS_HPP__
