#include <basic_marcher.hpp>
#include <olim.hpp>

#include <cmath>
#include <map>

#define REPS 3
#include "timer.hpp"

template <class marcher>
std::map<int, double>
time_marcher(int npowmin, int npowmax, speed_func s, int reps = REPS) {
  std::map<int, double> times;
  for (int npow = npowmin; npow <= npowmax; ++npow) {
    int n = (1 << npow) + 1, i = n/2;
    double h = 2./(n - 1);
    auto t = time([=] () {
      marcher m {n, n, h, s, 1., 1.};
      m.add_boundary_node(i, i);
      m.run();
    }, reps);
    times[n] = t;
  }
  return times;
}

enum class error_type {REL_L_INF, RMS};

template <class marcher>
double
rel_l_inf_error(int n, speed_func s, speed_func f) {
  int i0 = n/2;
  double h = 2./(n - 1);
  marcher m {n, n, h, s, 1., 1.};
  m.add_boundary_node(i0, i0);
  m.run();
  double abs_diff_max = 0;
  double u, u_max = 0;
  double U, U_max = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      u = f(h*j - 1, h*i - 1);
      U = m.get_value(i, j);
      u_max = std::max(u_max, u);
      U_max = std::max(U_max, U);
      abs_diff_max = std::max(fabs(u - U), abs_diff_max);
    }
  }
  return abs_diff_max/std::max(u_max, U_max);
}

template <class marcher>
double
rms_error(int n, speed_func s, speed_func f) {
  int i0 = n/2;
  double h = 2./(n - 1);
  marcher m {n, n, h, s, 1., 1.};
  m.add_boundary_node(i0, i0);
  m.run();
  double sum = 0, u, U;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      u = f(h*j - 1, h*i - 1);
      U = m.get_value(i, j);
      sum += (u - U)*(u - U);
    }
  }
  return std::sqrt(sum/std::pow(n, 2));
}

template <class marcher>
std::map<int, double>
marcher_errors(int npowmin, int npowmax, speed_func s, speed_func f,
               error_type error) {
  std::map<int, double> errors;
  for (int npow = npowmin; npow <= npowmax; ++npow) {
    int n = (1 << npow) + 1;
    if (error == error_type::REL_L_INF) {
      errors[n] = rel_l_inf_error<marcher>(n, s, f);
    } else {
      errors[n] = rms_error<marcher>(n, s, f);
    }
  }
  return errors;
}

#define GEN_CSV_DATA(MARCHER, S, F) do {                            \
    auto times = time_marcher<MARCHER>(npowmin, npowmax, S);        \
    auto l_inf_errors = marcher_errors<MARCHER>(                    \
      npowmin, npowmax, S, F, error_type::REL_L_INF);               \
    auto rms_errors = marcher_errors<MARCHER>(                      \
      npowmin, npowmax, S, F, error_type::RMS);                     \
    std::cout << #MARCHER << std::endl;                             \
    for (int npow = npowmin; npow <= npowmax; ++npow) {             \
      int n = (1 << npow) + 1;                                      \
      std::cout << n << ", " << times[n] << ", " << l_inf_errors[n] \
                << ", " << rms_errors[n] << std::endl;              \
    }                                                               \
  } while (0)

int main(int argc, char * argv[])
{
  speed_func s = s1, f = f1;

  if (argc >= 2) {
    std::string s_name = argv[1];
    if (s_name == "s0") {
      s = default_speed_func;
      f = default_speed_func_soln;
    } else if (s_name == "s1") {
      s = s1;
      f = f1;
    } else if (s_name == "s2") {
      s = s2;
      f = f2;
    } else if (s_name == "s3") {
      s = s3;
      f = f3;
    } else if (s_name == "s4xy") {
      // Doesn't work
      s = s4xy;
      f = f4xy;
    } else if (s_name == "s4xz") {
      // Doesn't work
      s = s4xz;
      f = f4xz;
    } else if (s_name == "s4yz") {
      // Doesn't work
      s = s4yz;
      f = f4yz;
    } else if (s_name == "s5") {
      s = s5;
      f = f5;
    } else if (s_name == "s6") {
      s = s6;
      f = f6;
    } else if (s_name == "s7") {
      s = s7;
      f = f7;
    }
  }

  int npowmin = 2, npowmax = 12;

  if (argc >= 3) npowmin = std::stoi(argv[2]);
  if (argc >= 4) npowmax = std::stoi(argv[3]);

  GEN_CSV_DATA(basic_marcher, s, f);
  GEN_CSV_DATA(olim4_rhr, s, f);
  GEN_CSV_DATA(olim4_mp0, s, f);
  GEN_CSV_DATA(olim4_mp1, s, f);
  GEN_CSV_DATA(olim8_rhr, s, f);
  GEN_CSV_DATA(olim8_mp0, s, f);
  GEN_CSV_DATA(olim8_mp1, s, f);
}
