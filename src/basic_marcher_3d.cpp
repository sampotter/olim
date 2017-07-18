#include "basic_marcher_3d.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#define COMPUTE_DISC() \
  (3*sh_sq - 2*(T1*T1 + T2*T2 + T3*T3 - T1*T2 - T1*T3 - T2*T3))

#define COMPUTE_VALUE() ((T1 + T2 + T3 + std::sqrt(disc))/3)

#define GET_VALUE(i) (nb[i]->get_value())

#define ADJ_TRI_UPDATE(dir1, dir2, dir3)                          \
  do {                                                            \
    if (has_octant[dir1 ## dir2 ## dir3]) {                       \
      T1 = GET_VALUE(DIR_ ## dir1), T2 = GET_VALUE(DIR_ ## dir2), \
        T3 = GET_VALUE(DIR_ ## dir3);                             \
      T = disc = COMPUTE_DISC() > 0 ? COMPUTE_VALUE() : T;        \
    }                                                             \
  } while (0);

enum neighbor {DIR_U, DIR_N, DIR_E, DIR_S, DIR_W, DIR_D};
enum octant {UNE, UES, USW, UWN, DNE, DES, DSW, DWN};

void basic_marcher_3d::update_node_value_impl(int i, int j, int k, double & T) {
  abstract_node * nb[6];
  get_valid_neighbors(i, j, k, nb);
  double sh = get_h()*S(i, j, k), sh_sq = sh*sh;
  double T1 = 0, T2 = 0, T3 = 0, disc = 0;

  bool has_nb[6] = {nb[0], nb[1], nb[2], nb[3], nb[4], nb[5]};
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

  // Triple point updates:

  // TODO: could replace with a macro using macro concatenation... (i.e. U ## N ## E)
  ADJ_TRI_UPDATE(U, N, E);
  ADJ_TRI_UPDATE(U, E, S);
  ADJ_TRI_UPDATE(U, S, W);
  ADJ_TRI_UPDATE(U, W, N);
  ADJ_TRI_UPDATE(D, N, E);
  ADJ_TRI_UPDATE(D, E, S);
  ADJ_TRI_UPDATE(D, S, W);
  ADJ_TRI_UPDATE(D, W, N);

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
