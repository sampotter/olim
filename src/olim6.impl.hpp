#ifndef __OLIM6_IMPL_HPP__
#define __OLIM6_IMPL_HPP__

#include <algorithm>

#include "common.macros.hpp"
#include "olim6.defs.hpp"

template <class update_rules>
void olim6<update_rules>::update_impl(int i, int j, int k, double & T) {
  abstract_node * nb[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  get_valid_neighbors(i, j, k, nb);
  double h = get_h(), s = S(i, j, k);
  
  double s_[6];
  for (int l = 0; l < 6; ++l) {
    if (nb[l]) {
      s_[l] = S(i + di[l], j + dj[l], k + dk[l]);
    }
  }

  bool quad[12] = {
    nb[DIR_N] && nb[DIR_E],
    nb[DIR_E] && nb[DIR_S],
    nb[DIR_S] && nb[DIR_W],
    nb[DIR_W] && nb[DIR_N],
    nb[DIR_U] && nb[DIR_N],
    nb[DIR_U] && nb[DIR_E],
    nb[DIR_U] && nb[DIR_S],
    nb[DIR_U] && nb[DIR_W],
    nb[DIR_D] && nb[DIR_N],
    nb[DIR_D] && nb[DIR_E],
    nb[DIR_D] && nb[DIR_S],
    nb[DIR_D] && nb[DIR_W]
  };

  bool oct[8] = {
    nb[DIR_U] && nb[DIR_N] && nb[DIR_E],
    nb[DIR_U] && nb[DIR_E] && nb[DIR_S],
    nb[DIR_U] && nb[DIR_S] && nb[DIR_W],
    nb[DIR_U] && nb[DIR_W] && nb[DIR_N],
    nb[DIR_D] && nb[DIR_N] && nb[DIR_E],
    nb[DIR_D] && nb[DIR_E] && nb[DIR_S],
    nb[DIR_D] && nb[DIR_S] && nb[DIR_W],
    nb[DIR_D] && nb[DIR_W] && nb[DIR_N]
  };

  // One point updates:

  // TODO: if we used a bit vector for the quadrants, we could just
  // use a mask to check this instead of using a larger boolean
  // expression---try this later, though... low priority

  if (!(quad[NE] || quad[ES] || quad[SW] || quad[WN])) {
    if (nb[DIR_U]) {
      T = std::min(T, this->adj1pt(GET_VALUE(DIR_U), s, s_[DIR_U], h));
    }
    if (nb[DIR_D]) {
      T = std::min(T, this->adj1pt(GET_VALUE(DIR_D), s, s_[DIR_D], h));
    }
  }
  if (!(quad[UE] || quad[DE] || quad[DW] || quad[UW])) {
    if (nb[DIR_N]) {
      T = std::min(T, this->adj1pt(GET_VALUE(DIR_N), s, s_[DIR_N], h));
    }
    if (nb[DIR_S]) {
      T = std::min(T, this->adj1pt(GET_VALUE(DIR_S), s, s_[DIR_S], h));
    }
  }
  if (!(quad[UN] || quad[DN] || quad[DS] || quad[US])) {
    if (nb[DIR_E]) {
      T = std::min(T, this->adj1pt(GET_VALUE(DIR_E), s, s_[DIR_E], h));
    }
    if (nb[DIR_W]) {
      T = std::min(T, this->adj1pt(GET_VALUE(DIR_W), s, s_[DIR_W], h));
    }
  }

  // Two point updates:

  if (!oct[UNE]) {
    if (!oct[UWN] && quad[UN]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_U), GET_VALUE(DIR_N), s, s_[DIR_U], s_[DIR_N], h));
    }
    if (!oct[DNE] && quad[NE]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_N), GET_VALUE(DIR_E), s, s_[DIR_N], s_[DIR_E], h));
    }
    if (!oct[UES] && quad[UE]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_U), GET_VALUE(DIR_E), s, s_[DIR_U], s_[DIR_E], h));
    }
  }
  if (!oct[USW]) {
    if (!oct[UES] && quad[US]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_U), GET_VALUE(DIR_S), s, s_[DIR_U], s_[DIR_S], h));
    }
    if (!oct[UWN] && quad[UW]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_U), GET_VALUE(DIR_W), s, s_[DIR_U], s_[DIR_W], h));
    }
    if (!oct[DSW] && quad[SW]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_S), GET_VALUE(DIR_W), s, s_[DIR_S], s_[DIR_W], h));
    }
  }
  if (!oct[DWN]) {
    if (!oct[UWN] && quad[WN]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_N), GET_VALUE(DIR_W), s, s_[DIR_N], s_[DIR_W], h));
    }
    if (!oct[DNE] && quad[DN]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_D), GET_VALUE(DIR_N), s, s_[DIR_D], s_[DIR_N], h));
    }
    if (!oct[DSW] && quad[DW]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_D), GET_VALUE(DIR_W), s, s_[DIR_D], s_[DIR_W], h));
    }
  }
  if (!oct[DES]) {
    if (!oct[DSW] && quad[DS]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_D), GET_VALUE(DIR_S), s, s_[DIR_D], s_[DIR_S], h));
    }
    if (!oct[DNE] && quad[DE]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_D), GET_VALUE(DIR_E), s, s_[DIR_D], s_[DIR_E], h));
    }
    if (!oct[UES] && quad[ES]) {
      T = std::min(
        T, this->adj2pt(
          GET_VALUE(DIR_E), GET_VALUE(DIR_S), s, s_[DIR_E], s_[DIR_S], h));
    }
  }

  // Three point updates:

  if (oct[UNE]) {
    T = std::min(
      T, this->adj3pt(
        GET_VALUE(DIR_U), GET_VALUE(DIR_N), GET_VALUE(DIR_E), s,
        s_[DIR_U], s_[DIR_N], s_[DIR_E], h));
  }
  if (oct[UES]) {
    T = std::min(
      T, this->adj3pt(
        GET_VALUE(DIR_U), GET_VALUE(DIR_E), GET_VALUE(DIR_S), s,
        s_[DIR_U], s_[DIR_E], s_[DIR_S], h));
  }
  if (oct[USW]) {
    T = std::min(
      T, this->adj3pt(
        GET_VALUE(DIR_U), GET_VALUE(DIR_S), GET_VALUE(DIR_W), s,
        s_[DIR_U], s_[DIR_S], s_[DIR_W], h));
  }
  if (oct[UWN]) {
    T = std::min(
      T, this->adj3pt(
        GET_VALUE(DIR_U), GET_VALUE(DIR_W), GET_VALUE(DIR_N), s,
        s_[DIR_U], s_[DIR_W], s_[DIR_N], h));
  }
  if (oct[DNE]) {
    T = std::min(
      T, this->adj3pt(
        GET_VALUE(DIR_D), GET_VALUE(DIR_N), GET_VALUE(DIR_E), s,
        s_[DIR_D], s_[DIR_N], s_[DIR_E], h));
  }
  if (oct[DES]) {
    T = std::min(
      T, this->adj3pt(
        GET_VALUE(DIR_D), GET_VALUE(DIR_E), GET_VALUE(DIR_S), s,
        s_[DIR_D], s_[DIR_E], s_[DIR_S], h));
  }
  if (oct[DSW]) {
    T = std::min(
      T, this->adj3pt(
        GET_VALUE(DIR_D), GET_VALUE(DIR_S), GET_VALUE(DIR_W), s,
        s_[DIR_D], s_[DIR_S], s_[DIR_W], h));
  }
  if (oct[DWN]) {
    T = std::min(
      T, this->adj3pt(
        GET_VALUE(DIR_D), GET_VALUE(DIR_W), GET_VALUE(DIR_N), s,
        s_[DIR_D], s_[DIR_W], s_[DIR_N], h));
  }
}

#endif // __OLIM6_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
