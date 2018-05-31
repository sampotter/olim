#include <basic_marcher_3d.hpp>
#include <olim.hpp>
#include <olim3d.hpp>

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <string>

#define REPS 5

using namespace std::chrono;

system_clock::time_point t0;

void tic() {
  t0 = system_clock::now();
}

double toc() {
  return duration<double, std::milli>(system_clock::now() - t0).count()/1000;
}

double time(std::function<void()> thunk, int reps = REPS) {
  double t = std::numeric_limits<double>::infinity();
  for (int i = 0; i < reps; ++i) {
    tic();
    thunk();
    t = std::min(t, toc());
  }
  return t;
}

template <class marcher_3d>
std::map<int, double>
time_marcher_3d(int npowmin, int npowmax, speed_func_3d s, int reps = REPS) {
  std::map<int, double> times;
  for (int npow = npowmin; npow <= npowmax; ++npow) {
    int n = (1 << npow) + 1, i = n/2;
    double h = 2./(n - 1);
    auto t = time([=] () {
      marcher_3d m {n, n, n, h, s, 1., 1., 1.};
      m.add_boundary_node(i, i, i);
      m.run();
    }, reps);
    times[n] = t;
  }
  return times;
}

enum class error_type {REL_L_INF, RMS};

template <class marcher_3d>
double
rel_l_inf_error(int n, speed_func_3d s, speed_func_3d f) {
  int i0 = n/2;
  double h = 2./(n - 1);
  marcher_3d m {n, n, n, h, s, 1., 1., 1.};
  m.add_boundary_node(i0, i0, i0);
  m.run();
  double abs_diff_max = 0;
  double u, u_max = 0;
  double U, U_max = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        u = f(h*j - 1, h*i - 1, h*k - 1);
        U = m.get_value(i, j, k);
        u_max = std::max(u_max, u);
        U_max = std::max(U_max, U);
        abs_diff_max = std::max(fabs(u - U), abs_diff_max);
      }
    }
  }
  return abs_diff_max/std::max(u_max, U_max);
}

template <class marcher_3d>
double
rms_error(int n, speed_func_3d s, speed_func_3d f) {
  int i0 = n/2;
  double h = 2./(n - 1);
  marcher_3d m {n, n, n, h, s, 1., 1., 1.};
  m.add_boundary_node(i0, i0, i0);
  m.run();
  double sum = 0, u, U;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        u = f(h*j - 1, h*i - 1, h*k - 1);
        U = m.get_value(i, j, k);
        sum += (u - U)*(u - U);
      }
    }
  }
  return std::sqrt(sum/std::pow(n, 3));
}

template <class marcher_3d>
std::map<int, double>
marcher_3d_errors(int npowmin, int npowmax, speed_func_3d s, speed_func_3d f,
                  error_type error) {
  std::map<int, double> errors;
  for (int npow = npowmin; npow <= npowmax; ++npow) {
    int n = (1 << npow) + 1;
    if (error == error_type::REL_L_INF) {
      errors[n] = rel_l_inf_error<marcher_3d>(n, s, f);
    } else {
      errors[n] = rms_error<marcher_3d>(n, s, f);
    }
  }
  return errors;
}

#define GEN_CSV_DATA(MARCHER, S3D, F3D) do {                        \
    auto times = time_marcher_3d<MARCHER>(npowmin, npowmax, S3D);   \
    auto l_inf_errors = marcher_3d_errors<MARCHER>(                 \
      npowmin, npowmax, S3D, F3D, error_type::REL_L_INF);           \
    auto rms_errors = marcher_3d_errors<MARCHER>(                   \
      npowmin, npowmax, S3D, F3D, error_type::RMS);                 \
    std::cout << #MARCHER << std::endl;                             \
    for (int npow = npowmin; npow <= npowmax; ++npow) {             \
      int n = (1 << npow) + 1;                                      \
      std::cout << n << ", " << times[n] << ", " << l_inf_errors[n] \
                << ", " << rms_errors[n] << std::endl;              \
    }                                                               \
  } while (0)

int main(int argc, char * argv[]) {
  (void) argc;
  (void) argv;

  int npowmin = 2, npowmax = 6;

  // std::string speed_func_names[] = {
  //   "s0", "s1", "s2", "s3", "s4", "s5", "s6", "s7"
  // };

  // speed_func_3d s3ds[] = {
  //   default_speed_func, s1, s2, s3, s4, s5, s6, s7
  // };

  // speed_func_3d f3ds[] = {
  //   default_speed_func_soln, f1, f2, f3, f4, f5, f6, f7
  // };

  // speed_func_3d s3d = default_speed_func, f3d = default_speed_func_soln;
  // speed_func_3d s3d = s4, f3d = f4;

  // for (int i = 0; i < sizeof(s3ds)/sizeof(speed_func_3d); ++i) {
  //   speed_func_3d s3d = s3ds[i], f3d = f3ds[i];

  //   std::cout << speed_func_names[i] << std::endl;

    GEN_CSV_DATA(basic_marcher_3d, s3d, f3d);
    GEN_CSV_DATA(olim6_mp0, s3d, f3d);
    GEN_CSV_DATA(olim6_mp1, s3d, f3d);
    GEN_CSV_DATA(olim6_rhr, s3d, f3d);
    GEN_CSV_DATA(olim18_mp0, s3d, f3d);
    GEN_CSV_DATA(olim18_mp1, s3d, f3d);
    GEN_CSV_DATA(olim18_rhr, s3d, f3d);
    GEN_CSV_DATA(olim26_mp0, s3d, f3d);
    GEN_CSV_DATA(olim26_mp1, s3d, f3d);
    GEN_CSV_DATA(olim26_rhr, s3d, f3d);
    GEN_CSV_DATA(olim3d_hu_mp0, s3d, f3d);
    GEN_CSV_DATA(olim3d_hu_mp1, s3d, f3d);
    GEN_CSV_DATA(olim3d_hu_rhr, s3d, f3d);
// }
}
