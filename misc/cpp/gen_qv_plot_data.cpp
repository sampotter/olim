#include <basic_marcher.hpp>
#include <basic_marcher_3d.hpp>
#include <olim.hpp>
#include <olim3d.hpp>

#include <cstdio>
#include <string>

#define FACTORED 1

constexpr double a = 2, vx = 5, vy = 20, vz = 13;
constexpr double X1 = 0, Y1 = 0.5, Z1 = 0;
constexpr double X2 = 0.8, Y2 = 0.2, Z2 = 0;
constexpr double X3 = 0.6, Y3 = 1, Z3 = 0;
constexpr double r = 0.01;

constexpr int minpow2d = 3;
constexpr int maxpow2d = 7;

constexpr int minpow3d = 3;
constexpr int maxpow3d = 6;

const double v2rn = 1/sqrt(vx*vx + vy*vy);
const double v3rn = 1/sqrt(vx*vx + vy*vy + vz*vz);

double qv_s_2d(double x, double y) {
  return 1./(a + vx*x + vy*y);
}

double qv_u_2d(double x, double y) {
  double s1 = qv_s_2d(X1, Y1);
  double u1 = v2rn*acosh(
    1.0 + s1*qv_s_2d(x, y)*(vx*vx + vy*vy)*(
      (x - X1)*(x - X1) + (y - Y1)*(y - Y1))/2);

  double s2 = qv_s_2d(X2, Y2);
  double u2 = v2rn*acosh(
    1.0 + s2*qv_s_2d(x, y)*(vx*vx + vy*vy)*(
      (x - X2)*(x - X2) + (y - Y2)*(y - Y2))/2);

  return std::min(u1, u2);
}

template <class olim>
void add_children_2d(olim & o, int i0, int j0, int n, double h) {
  for (int i = 0; i < n; ++i) {
    double y = h*i;
    for (int j = 0; j < n; ++j) {
      double x = h*j;
      if (hypot(x, y) <= r) {
        o.set_node_fac_parent(i, j, i0, j0);
      }
    }
  }
}

template <class olim>
std::pair<double, double> do_qv_2d(int n) {
  int i1 = 0, j1 = n*X1;
  int i2 = 0, j2 = n*X2;
  double h = 1./(n - 1.);

  olim o {n, n, h, qv_s_2d};
  o.add_boundary_node(i1, j1);
  o.add_boundary_node(i2, j2);
#if FACTORED
  add_children_2d(o, i1, j1, n, h);
  add_children_2d(o, i2, j2, n, h);
#endif
  o.run();

  double rms, linf, u, U;

  // RMS error
  double sum = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      u = qv_u_2d(h*j, h*i);
      U = o.get_value(i, j);
      sum += (u - U)*(u - U);
    }
  }
  rms = sqrt(sum)/n;

  // rel linf error
  double umax = 0, Umax = 0, dmax = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      u = qv_u_2d(h*j, h*i);
      U = o.get_value(i, j);
      umax = std::max(umax, u);
      Umax = std::max(Umax, U);
      dmax = std::max(dmax, fabs(u - U));
    }
  }
  linf = dmax/std::max(umax, Umax);

  return {rms, linf};
}

template <class olim>
void print_errors_2d(std::string const & oname) {
  printf("%s\n", oname.c_str());
  double rms, linf;
  for (int p = minpow2d; p <= maxpow2d; ++p) {
    int n = (1 << p) + 1;
    std::tie(rms, linf) = do_qv_2d<olim>(n);
    printf("%d, %g, %g\n", n, rms, linf);
  }
}

double qv_s_3d(double x, double y, double z) {
  return 1./(a + vx*x + vy*y + vz*z);
}

double qv_u_3d(double x, double y, double z) {
  double s1 = qv_s_3d(X1, Y1, Y1); // TODO: double check
  double u1 = v3rn*acosh(
    1.0 + s1*qv_s_3d(x, y, z)*(vx*vx + vy*vy + vz*vz)*(
      (x - X1)*(x - X1) + (y - Y1)*(y - Y1) + (z - Z1)*(z - X1))/2);

  double s2 = qv_s_3d(X2, Y2, Y2); // TODO: double check
  double u2 =  v3rn*acosh(
    1.0 + s2*qv_s_3d(x, y, z)*(vx*vx + vy*vy + vz*vz)*(
      (x - X2)*(x - X2) + (y - Y2)*(y - Y2) + (z - Z2)*(z - X2))/2);

  double s3 = qv_s_3d(X3, Y3, Y3); // TODO: double check
  double u3 = v3rn*acosh(
    1.0 + s3*qv_s_3d(x, y, z)*(vx*vx + vy*vy + vz*vz)*(
      (x - X3)*(x - X3) + (y - Y3)*(y - Y3) + (z - Z3)*(z - X3))/2);

  return std::min(u1, std::min(u2, u3));
}

template <class olim>
void add_children_3d(olim & o, int i0, int j0, int k0, double h, int n) {
  for (int i = 0; i < n; ++i) {
    double y = h*i;
    for (int j = 0; j < n; ++j) {
      double x = h*j;
      for (int k = 0; k < n; ++k) {
        double z = h*k;
        if (sqrt(x*x + y*y* + z*z) <= r) {
          o.set_node_fac_parent(i, j, k, i0, j0, k0);
        }
      }
    }
  }
}

template <class olim>
std::pair<double, double> do_qv_3d(int n) {
  int i1 = n*Y1, j1 = n*X1, k1 = 0;
  int i2 = n*Y2, j2 = n*X2, k2 = 0;
  int i3 = n*Y3, j3 = n*X3, k3 = 0;
  double h = 1./(n - 1.);

  olim o {n, n, n, h, qv_s_3d};
  o.add_boundary_node(i1, j1, k1);
  o.add_boundary_node(i2, j2, k2);
  o.add_boundary_node(i3, j3, k3);
#if FACTORED
  add_children_3d(o, i1, j1, k1, n, h);
  add_children_3d(o, i2, j2, k2, n, h);
  add_children_3d(o, i3, j3, k3, n, h);
#endif
  o.run();

  double rms, linf, u, U;

  // RMS error
  double sum = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        u = qv_u_3d(h*j, h*i, h*k);
        U = o.get_value(i, j, k);
        sum += (u - U)*(u - U);
      }
    }
  }
  rms = sqrt(sum/pow(n, 3));

  // rel linf error
  double umax = 0, Umax = 0, dmax = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        u = qv_u_3d(h*j, h*i, h*k);
        U = o.get_value(i, j, k);
        umax = std::max(umax, u);
        Umax = std::max(Umax, U);
        dmax = std::max(dmax, fabs(u - U));
      }
    }
  }
  linf = dmax/std::max(umax, Umax);

  return {rms, linf};
}

template <class olim>
void print_errors_3d(std::string const & oname) {
  printf("%s\n", oname.c_str());
  double rms, linf;
  for (int p = minpow3d; p <= maxpow3d; ++p) {
    int n = (1 << p) + 1;
    std::tie(rms, linf) = do_qv_3d<olim>(n);
    printf("%d, %g, %g\n", n, rms, linf);
  }
}

// int main(int argc, char * argv[]) {
int main() {
  print_errors_2d<olim4_rhr>("olim4_rhr");
  print_errors_2d<olim4_mp0>("olim4_mp0");
  print_errors_2d<olim4_mp1>("olim4_mp1");
  print_errors_2d<olim8_rhr>("olim8_rhr");
  print_errors_2d<olim8_mp0>("olim8_mp0");
  print_errors_2d<olim8_mp1>("olim8_mp1");
  print_errors_3d<basic_marcher_3d>("basic_marcher_3d");
  print_errors_3d<olim6_rhr>("olim6_rhr");
  print_errors_3d<olim6_mp0>("olim6_mp0");
  print_errors_3d<olim6_mp1>("olim6_mp1");
  print_errors_3d<olim18_rhr>("olim18_rhr");
  print_errors_3d<olim18_mp0>("olim18_mp0");
  print_errors_3d<olim18_mp1>("olim18_mp1");
  print_errors_3d<olim26_rhr>("olim26_rhr");
  print_errors_3d<olim26_mp0>("olim26_mp0");
  print_errors_3d<olim26_mp1>("olim26_mp1");
  print_errors_3d<olim3d_hu_rhr>("olim3d_hu_rhr");
  print_errors_3d<olim3d_hu_mp0>("olim3d_hu_mp0");
  print_errors_3d<olim3d_hu_mp1>("olim3d_hu_mp1");
}
