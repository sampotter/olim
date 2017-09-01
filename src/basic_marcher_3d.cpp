#include "basic_marcher_3d.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "common.macros.hpp"
#include "olim6.defs.hpp"

#define COMPUTE_DISC_2PT() (2*sh*sh - (T1 - T2)*(T1 - T2))

#define COMPUTE_VALUE_2PT() ((T1 + T2 + std::sqrt(disc))/2)

#define COMPUTE_DISC_3PT() \
  (3*sh_sq - 2*(T1*T1 + T2*T2 + T3*T3 - T1*T2 - T1*T3 - T2*T3))

#define COMPUTE_VALUE_3PT() ((T1 + T2 + T3 + std::sqrt(disc))/3)

void basic_marcher_3d::update_impl(int i, int j, int k, double & T) {
  abstract_node * nb[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  get_valid_neighbors(i, j, k, nb);
  double sh = get_h()*S(i, j, k), sh_sq = sh*sh;
  double T1 = 0, T2 = 0, T3 = 0, disc = 0;

  bool has_nb[6] = {nb[0], nb[1], nb[2], nb[3], nb[4], nb[5]};

  bool has_quadrant[12] = {
    has_nb[DIR_N] && has_nb[DIR_E],
    has_nb[DIR_E] && has_nb[DIR_S],
    has_nb[DIR_S] && has_nb[DIR_W],
    has_nb[DIR_W] && has_nb[DIR_N],
    has_nb[DIR_U] && has_nb[DIR_N],
    has_nb[DIR_U] && has_nb[DIR_E],
    has_nb[DIR_U] && has_nb[DIR_S],
    has_nb[DIR_U] && has_nb[DIR_W],
    has_nb[DIR_D] && has_nb[DIR_N],
    has_nb[DIR_D] && has_nb[DIR_E],
    has_nb[DIR_D] && has_nb[DIR_S],
    has_nb[DIR_D] && has_nb[DIR_W]
  };

  bool has_octant[8] = {
    has_nb[DIR_U] && has_nb[DIR_N] && has_nb[DIR_E],
    has_nb[DIR_U] && has_nb[DIR_E] && has_nb[DIR_S],
    has_nb[DIR_U] && has_nb[DIR_S] && has_nb[DIR_W],
    has_nb[DIR_U] && has_nb[DIR_W] && has_nb[DIR_N],
    has_nb[DIR_D] && has_nb[DIR_N] && has_nb[DIR_E],
    has_nb[DIR_D] && has_nb[DIR_E] && has_nb[DIR_S],
    has_nb[DIR_D] && has_nb[DIR_S] && has_nb[DIR_W],
    has_nb[DIR_D] && has_nb[DIR_W] && has_nb[DIR_N]
  };

  // Two point updates:

  if (has_quadrant[NE]) {
    T1 = GET_VALUE(DIR_N), T2 = GET_VALUE(DIR_E);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[ES]) {
    T1 = GET_VALUE(DIR_E), T2 = GET_VALUE(DIR_S);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[SW]) {
    T1 = GET_VALUE(DIR_S), T2 = GET_VALUE(DIR_W);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[WN]) {
    T1 = GET_VALUE(DIR_W), T2 = GET_VALUE(DIR_N);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[UN]) {
    T1 = GET_VALUE(DIR_U), T2 = GET_VALUE(DIR_N);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[UE]) {
    T1 = GET_VALUE(DIR_U), T2 = GET_VALUE(DIR_E);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[US]) {
    T1 = GET_VALUE(DIR_U), T2 = GET_VALUE(DIR_S);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[UW]) {
    T1 = GET_VALUE(DIR_U), T2 = GET_VALUE(DIR_W);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[DN]) {
    T1 = GET_VALUE(DIR_D), T2 = GET_VALUE(DIR_N);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[DE]) {
    T1 = GET_VALUE(DIR_D), T2 = GET_VALUE(DIR_E);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[DS]) {
    T1 = GET_VALUE(DIR_D), T2 = GET_VALUE(DIR_S);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (has_quadrant[DW]) {
    T1 = GET_VALUE(DIR_D), T2 = GET_VALUE(DIR_W);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_2PT());
    }
  }

  // Triple point updates:

  if (has_octant[UNE]) {
    T1 = GET_VALUE(DIR_U), T2 = GET_VALUE(DIR_N), T3 = GET_VALUE(DIR_E);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (has_octant[UES]) {
    T1 = GET_VALUE(DIR_U), T2 = GET_VALUE(DIR_E), T3 = GET_VALUE(DIR_S);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (has_octant[USW]) {
    T1 = GET_VALUE(DIR_U), T2 = GET_VALUE(DIR_S), T3 = GET_VALUE(DIR_W);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (has_octant[UWN]) {
    T1 = GET_VALUE(DIR_U), T2 = GET_VALUE(DIR_W), T3 = GET_VALUE(DIR_N);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (has_octant[DNE]) {
    T1 = GET_VALUE(DIR_D), T2 = GET_VALUE(DIR_N), T3 = GET_VALUE(DIR_E);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (has_octant[DES]) {
    T1 = GET_VALUE(DIR_D), T2 = GET_VALUE(DIR_E), T3 = GET_VALUE(DIR_S);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (has_octant[DSW]) {
    T1 = GET_VALUE(DIR_D), T2 = GET_VALUE(DIR_S), T3 = GET_VALUE(DIR_W);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (has_octant[DWN]) {
    T1 = GET_VALUE(DIR_D), T2 = GET_VALUE(DIR_W), T3 = GET_VALUE(DIR_N);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = std::min(T, COMPUTE_VALUE_3PT());
    }
  }

  // Single point updates:

  // TODO: could order this so that fewer comparisons are made
  if (has_nb[DIR_U] && !(has_octant[UNE] || has_octant[UES] ||
                         has_octant[USW] || has_octant[UWN])) {
    T = std::min(T, GET_VALUE(DIR_U) + sh);
  }
  if (has_nb[DIR_N] && !(has_octant[UNE] || has_octant[UWN] ||
                         has_octant[DNE] || has_octant[DWN])) {
    T = std::min(T, GET_VALUE(DIR_N) + sh);
  }
  if (has_nb[DIR_E] && !(has_octant[UNE] || has_octant[UES] ||
                         has_octant[DNE] || has_octant[DES])) {
    T = std::min(T, GET_VALUE(DIR_E) + sh);
  }
  if (has_nb[DIR_S] && !(has_octant[UES] || has_octant[USW] ||
                         has_octant[DES] || has_octant[DSW])) {
    T = std::min(T, GET_VALUE(DIR_S) + sh);
  }
  if (has_nb[DIR_W] && !(has_octant[USW] || has_octant[UWN] ||
                         has_octant[DSW] || has_octant[DWN])) {
    T = std::min(T, GET_VALUE(DIR_W) + sh);
  }
  if (has_nb[DIR_D] && !(has_octant[DNE] || has_octant[DWN] ||
                         has_octant[DSW] || has_octant[DWN])) {
    T = std::min(T, GET_VALUE(DIR_D) + sh);
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
