#ifndef __OLIM6_IMPL_HPP__
#define __OLIM6_IMPL_HPP__

#include <algorithm>

#include "common.macros.hpp"
#include "olim6.defs.hpp"

template <class update_rules>
void olim6<update_rules>::update_impl(int i, int j, int k, double & T) {
  using namespace olim6_defs;
  using std::min;

  abstract_node * nb[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  get_valid_neighbors(i, j, k, nb);
  double h = get_h(), s = speed(i, j, k);
  
  double s_[6];
  for (int l = 0; l < 6; ++l) {
    if (nb[l]) {
      s_[l] = speed(i + di[l], j + dj[l], k + dk[l]);
    }
  }

  bool quad[12] = {
    nb[N] && nb[E],
    nb[E] && nb[S],
    nb[S] && nb[W],
    nb[W] && nb[N],
    nb[U] && nb[N],
    nb[U] && nb[E],
    nb[U] && nb[S],
    nb[U] && nb[W],
    nb[D] && nb[N],
    nb[D] && nb[E],
    nb[D] && nb[S],
    nb[D] && nb[W]
  };

  bool oct[8] = {
    nb[U] && nb[N] && nb[E],
    nb[U] && nb[E] && nb[S],
    nb[U] && nb[S] && nb[W],
    nb[U] && nb[W] && nb[N],
    nb[D] && nb[N] && nb[E],
    nb[D] && nb[E] && nb[S],
    nb[D] && nb[S] && nb[W],
    nb[D] && nb[W] && nb[N]
  };

  // One point updates:

  // TODO: if we used a bit vector for the quadrants, we could just
  // use a mask to check this instead of using a larger boolean
  // expression---try this later, though... low priority

  if (!(quad[NE] || quad[ES] || quad[SW] || quad[WN])) {
    if (nb[U]) T = min(T, this->adj1pt(VAL(U), s, s_[U], h));
    if (nb[D]) T = min(T, this->adj1pt(VAL(D), s, s_[D], h));
  }
  if (!(quad[UE] || quad[DE] || quad[DW] || quad[UW])) {
    if (nb[N]) T = min(T, this->adj1pt(VAL(N), s, s_[N], h));
    if (nb[S]) T = min(T, this->adj1pt(VAL(S), s, s_[S], h));
  }
  if (!(quad[UN] || quad[DN] || quad[DS] || quad[US])) {
    if (nb[E]) T = min(T, this->adj1pt(VAL(E), s, s_[E], h));
    if (nb[W]) T = min(T, this->adj1pt(VAL(W), s, s_[W], h));
  }

  // Two point updates:

  if (!oct[UNE]) {
    if (!oct[UWN] && quad[UN]) {
      T = min(T, this->adj2pt(VAL(U), VAL(N), s, s_[U], s_[N], h));
    }
    if (!oct[DNE] && quad[NE]) {
      T = min(T, this->adj2pt(VAL(N), VAL(E), s, s_[N], s_[E], h));
    }
    if (!oct[UES] && quad[UE]) {
      T = min(T, this->adj2pt(VAL(U), VAL(E), s, s_[U], s_[E], h));
    }
  }
  if (!oct[USW]) {
    if (!oct[UES] && quad[US]) {
      T = min(T, this->adj2pt(VAL(U), VAL(S), s, s_[U], s_[S], h));
    }
    if (!oct[UWN] && quad[UW]) {
      T = min(T, this->adj2pt(VAL(U), VAL(W), s, s_[U], s_[W], h));
    }
    if (!oct[DSW] && quad[SW]) {
      T = min(T, this->adj2pt(VAL(S), VAL(W), s, s_[S], s_[W], h));
    }
  }
  if (!oct[DWN]) {
    if (!oct[UWN] && quad[WN]) {
      T = min(T, this->adj2pt(VAL(N), VAL(W), s, s_[N], s_[W], h));
    }
    if (!oct[DNE] && quad[DN]) {
      T = min(T, this->adj2pt(VAL(D), VAL(N), s, s_[D], s_[N], h));
    }
    if (!oct[DSW] && quad[DW]) {
      T = min(T, this->adj2pt(VAL(D), VAL(W), s, s_[D], s_[W], h));
    }
  }
  if (!oct[DES]) {
    if (!oct[DSW] && quad[DS]) {
      T = min(T, this->adj2pt(VAL(D), VAL(S), s, s_[D], s_[S], h));
    }
    if (!oct[DNE] && quad[DE]) {
      T = min(T, this->adj2pt(VAL(D), VAL(E), s, s_[D], s_[E], h));
    }
    if (!oct[UES] && quad[ES]) {
      T = min(T, this->adj2pt(VAL(E), VAL(S), s, s_[E], s_[S], h));
    }
  }

  // Three point updates:

  if (oct[UNE]) {
    T = min(T, this->adj3pt(VAL(U), VAL(N), VAL(E), s, s_[U], s_[N], s_[E], h));
  }
  if (oct[UES]) {
    T = min(T, this->adj3pt(VAL(U), VAL(E), VAL(S), s, s_[U], s_[E], s_[S], h));
  }
  if (oct[USW]) {
    T = min(T, this->adj3pt(VAL(U), VAL(S), VAL(W), s, s_[U], s_[S], s_[W], h));
  }
  if (oct[UWN]) {
    T = min(T, this->adj3pt(VAL(U), VAL(W), VAL(N), s, s_[U], s_[W], s_[N], h));
  }
  if (oct[DNE]) {
    T = min(T, this->adj3pt(VAL(D), VAL(N), VAL(E), s, s_[D], s_[N], s_[E], h));
  }
  if (oct[DES]) {
    T = min(T, this->adj3pt(VAL(D), VAL(E), VAL(S), s, s_[D], s_[E], s_[S], h));
  }
  if (oct[DSW]) {
    T = min(T, this->adj3pt(VAL(D), VAL(S), VAL(W), s, s_[D], s_[S], s_[W], h));
  }
  if (oct[DWN]) {
    T = min(T, this->adj3pt(VAL(D), VAL(W), VAL(N), s, s_[D], s_[W], s_[N], h));
  }
}

#endif // __OLIM6_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
