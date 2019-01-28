#pragma once

#include <array>

struct no_slow_t {};
constexpr no_slow_t no_slow {};

double s0(double x, double y);
double f0(double x, double y);

double s0(double x, double y, double z);
double f0(double x, double y, double z);

double s1(double x, double y);
double f1(double x, double y);

double s1(double x, double y, double z);
double f1(double x, double y, double z);

double s2(double x, double y);
double f2(double x, double y);

double s2(double x, double y, double z);
double f2(double x, double y, double z);

double s5(double x, double y);
double f5(double x, double y);

// Masha's slowness function
double s7(double x, double y);
double f7(double x, double y);

#ifndef __clang__
#    pragma GCC diagnostic ignored "-Wunused-variable"
#endif

static std::array<double (*)(double, double), 5> slow2s {{
  s0, s1, s2, s5, s7
}};

static std::array<double (*)(double, double), 5> soln2s {{
  f0, f1, f2, f5, f7
}};

static std::array<double (*)(double, double, double), 3> slow3s {{
  s0, s1, s2
}};

static std::array<double (*)(double, double, double), 3> soln3s {{
  f0, f1, s2
}};

#ifndef __clang__
#    pragma GCC diagnostic pop
#endif
