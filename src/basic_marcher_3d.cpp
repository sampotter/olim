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
  using namespace olim6_defs;
  using std::min;

  abstract_node * nb[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  get_valid_neighbors(i, j, k, nb);
  double sh = get_h()*speed(i, j, k), sh_sq = sh*sh;
  double T1 = 0, T2 = 0, T3 = 0, disc = 0;

  bool quad[12] = {
    nb[N] && nb[E], nb[E] && nb[S], nb[S] && nb[W], nb[W] && nb[N],
    nb[U] && nb[N], nb[U] && nb[E], nb[U] && nb[S], nb[U] && nb[W],
    nb[D] && nb[N], nb[D] && nb[E], nb[D] && nb[S], nb[D] && nb[W]
  };

  bool oct[8] = {
    nb[U] && nb[N] && nb[E], nb[U] && nb[E] && nb[S],
    nb[U] && nb[S] && nb[W], nb[U] && nb[W] && nb[N],
    nb[D] && nb[N] && nb[E], nb[D] && nb[E] && nb[S],
    nb[D] && nb[S] && nb[W], nb[D] && nb[W] && nb[N]
  };

  // Two point updates:

  if (quad[NE]) {
    T1 = VAL(N), T2 = VAL(E);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[ES]) {
    T1 = VAL(E), T2 = VAL(S);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[SW]) {
    T1 = VAL(S), T2 = VAL(W);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[WN]) {
    T1 = VAL(W), T2 = VAL(N);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[UN]) {
    T1 = VAL(U), T2 = VAL(N);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[UE]) {
    T1 = VAL(U), T2 = VAL(E);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[US]) {
    T1 = VAL(U), T2 = VAL(S);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[UW]) {
    T1 = VAL(U), T2 = VAL(W);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[DN]) {
    T1 = VAL(D), T2 = VAL(N);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[DE]) {
    T1 = VAL(D), T2 = VAL(E);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[DS]) {
    T1 = VAL(D), T2 = VAL(S);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }
  if (quad[DW]) {
    T1 = VAL(D), T2 = VAL(W);
    disc = COMPUTE_DISC_2PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_2PT());
    }
  }

  // Triple point updates:

  if (oct[UNE]) {
    T1 = VAL(U), T2 = VAL(N), T3 = VAL(E);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (oct[UES]) {
    T1 = VAL(U), T2 = VAL(E), T3 = VAL(S);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (oct[USW]) {
    T1 = VAL(U), T2 = VAL(S), T3 = VAL(W);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (oct[UWN]) {
    T1 = VAL(U), T2 = VAL(W), T3 = VAL(N);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (oct[DNE]) {
    T1 = VAL(D), T2 = VAL(N), T3 = VAL(E);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (oct[DES]) {
    T1 = VAL(D), T2 = VAL(E), T3 = VAL(S);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (oct[DSW]) {
    T1 = VAL(D), T2 = VAL(S), T3 = VAL(W);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_3PT());
    }
  }
  if (oct[DWN]) {
    T1 = VAL(D), T2 = VAL(W), T3 = VAL(N);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) {
      T = min(T, COMPUTE_VALUE_3PT());
    }
  }

  // Single point updates:

  // TODO: could order this so that fewer comparisons are made
  if (nb[U] && !(oct[UNE] || oct[UES] || oct[USW] || oct[UWN])) {
    T = min(T, VAL(U) + sh);
  }
  if (nb[N] && !(oct[UNE] || oct[UWN] || oct[DNE] || oct[DWN])) {
    T = min(T, VAL(N) + sh);
  }
  if (nb[E] && !(oct[UNE] || oct[UES] || oct[DNE] || oct[DES])) {
    T = min(T, VAL(E) + sh);
  }
  if (nb[S] && !(oct[UES] || oct[USW] || oct[DES] || oct[DSW])) {
    T = min(T, VAL(S) + sh);
  }
  if (nb[W] && !(oct[USW] || oct[UWN] || oct[DSW] || oct[DWN])) {
    T = min(T, VAL(W) + sh);
  }
  if (nb[D] && !(oct[DNE] || oct[DWN] || oct[DSW] || oct[DWN])) {
    T = min(T, VAL(D) + sh);
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
