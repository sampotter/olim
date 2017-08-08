#include "rootfinding.hpp"

#include <cmath>

rootfinder::rootfinder(double (* f)(double)): _f {f} {}

double rootfinder::eval_f(double x) const {
  return _f(x);
}

double secant::find_root(double x0, double x1, double tol) const {
  double x[3], f[3], dx;
  x[1] = x0;
  x[2] = x1;
  f[1] = eval_f(x[1]);
  f[2] = eval_f(x[2]);
  do {
    x[0] = (x[2]*f[1] - x[1]*f[2])/(f[1] - f[2]);
    dx = x[0] - x[1];
    x[2] = x[1];
    x[1] = x[0];
    f[2] = f[1];
    f[1] = eval_f(x[1]);
  } while (std::fabs(dx) > tol);
  return x[0];
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
