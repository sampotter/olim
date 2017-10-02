#ifndef __OLIM26_RECT_IMPL_HPP__
#define __OLIM26_RECT_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.macros.hpp"
#include "olim26.defs.hpp"

// neighbor order:
// degree 1: N, E, U, S, W, D
//           0, 1, 2, 3, 4, 5
// degree 2: UN, UE, US, UW, NE, SE, SW, NW, DN, DE, DS, DW
//           6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17
// degree 3: UNE, USE, USW, UNW, DNE, DSE, DSW, DNW
//           18,  19,  20,  21,  22,  23,  24,  25

template <class node, class update_rules>
int olim26_rect<node, update_rules>::di[] = {
  1, 0, 0, -1, 0, 0,
  1, 0, -1, 0, 1, -1, -1, 1, 1, 0, -1, 0,
  1, -1, -1, 1, 1, -1, -1, 1
};

template <class node, class update_rules>
int olim26_rect<node, update_rules>::dj[] = {
  0, 1, 0, 0, -1, 0,
  0, 1, 0, -1, 1, 1, -1, -1, 0, 1, 0, -1,
  1, 1, -1, -1, 1, 1, -1, -1
};

template <class node, class update_rules>
int olim26_rect<node, update_rules>::dk[] = {
  0, 0, 1, 0, 0, -1,
  1, 1, 1, 1, 0, 0, 0, 0, -1, -1, -1, -1,
  1, 1, 1, 1, -1, -1, -1, -1
};

template <class node, class update_rules>
void olim26_rect<node, update_rules>::get_valid_neighbors(int i, int j, int k,
                                                          abstract_node ** nb) {
  int a, b, c;
  for (int l = 0; l < 26; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->is_valid(a, b, c)) {
      nb[l] = &this->operator()(a, b, c);
    }
  }
}

template <class node, class update_rules>
void olim26_rect<node, update_rules>::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();
  int k = static_cast<node *>(n)->get_k();

#if PRINT_UPDATES
  printf("olim26_rect::stage_neighbors_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  for (int l = 0; l < 26; ++l) {
    this->stage(i + di[l], j + dj[l], k + dk[l]);
  }

  int a, b, c;
  for (int l = 0; l < 26; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      this->update(a, b, c);
    }
  }
}

namespace olim26 {
  int line1tris[6][8] = {
    {UN, UNE, NE, DNE, DN, DNW, NW, UNW}, // N
    {UE, USE, SE, DSE, DE, DNE, NE, UNE}, // E
    {UN, UNE, UE, USE, US, USW, UW, UNW}, // U
    {US, USE, SE, DSE, DS, DSW, SW, USW}, // S
    {UW, USW, SW, DSW, DW, DNW, NW, UNW}, // W
    {DN, DNE, DE, DSE, DS, DSW, DW, DNW}, // D
  };

  int line2tris[12][4] = {
    {U, N, UNE, UNW}, // UN
    {U, E, UNE, USE}, // UE
    {U, S, USE, USW}, // US
    {U, W, UNW, USW}, // UW
    {N, E, UNE, DNE}, // NE
    {S, E, USE, DSE}, // SE
    {S, W, USW, DSW}, // SW
    {N, W, UNW, DNW}, // NW
    {D, N, DNE, DNW}, // DN
    {D, E, DNE, DSE}, // DE
    {D, S, DSE, DSW}, // DS
    {D, W, DNW, DSW}, // DW
  };

  int line3tris[8][6] = {
    {U, UN, N, NE, E, UE}, // UNE
    {U, US, S, SE, E, UE}, // USE
    {U, UW, W, SW, S, US}, // USW
    {U, UN, N, NW, W, UW}, // UNW
    {D, DN, N, NE, E, DE}, // DNE
    {D, DS, S, SE, E, DE}, // DSE
    {D, DS, S, SW, W, DW}, // DSW
    {D, DN, N, NW, W, DW}, // DNW
  };
}

template <class node, class update_rules>
void olim26_rect<node, update_rules>::update_impl(int i, int j, int k, double & T) {
  using namespace olim26;
  using std::min;

#ifdef PRINT_UPDATES
  printf("olim26_rect::update_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  abstract_node * nb[26];
  memset(nb, 0x0, 26*sizeof(abstract_node *));
  get_valid_neighbors(i, j, k, nb);

  double h = this->get_h(), s = this->speed(i, j, k);
  int l, l0, l1, l2, a, b, * is;

  /**
   * Degree 1 line updates
   */
  for (l = 0; l < 6; ++l) {
    is = olim26::line1tris[l];
    if (!nb[l])
      continue;
    if (nb[is[0]] || nb[is[1]] || nb[is[2]] || nb[is[3]] ||
        nb[is[4]] || nb[is[5]] || nb[is[6]] || nb[is[7]])
      continue;
    T = min(T, this->line1(VAL(l), s, h));
  }

  /**
   * Degree 2 line updates
   */
  for (; l < 18; ++l) {
    is = olim26::line2tris[l - 6];
    if (!nb[l])
      continue;
    if (nb[is[0]] || nb[is[1]] || nb[is[2]] || nb[is[3]])
      continue;
    T = min(T, this->line2(VAL(l), s, h));
  }

  /**
   * Degree 3 line updates
   */
  for (; l < 26; ++l) {
    is = olim26::line3tris[l - 18];
    if (!nb[l])
      continue;
    if (nb[is[0]] || nb[is[1]] || nb[is[2]] ||
        nb[is[3]] || nb[is[4]] || nb[is[5]])
      continue;
    T = min(T, this->line3(VAL(l), s, h));
  }

  /**
   * Degree (1, 2) and (1, 3) triangle updates
   */
  for (l0 = 0; l0 < 6; ++l0) {
    is = olim26::line1tris[l0];
    if (nb[l0]) {
      for (a = 0, l1 = *is; a < 8; l1 = is[++a]) {
        if (nb[l1]) {
          if (a % 2 == 0) {
            T = min(T, this->tri12(VAL(l0), VAL(l1), s, h));
          } else {
            T = min(T, this->tri13(VAL(l1), VAL(l0), s, h));
          }
        }
      }
    }
  }

  /**
   * Degree (2, 3) triangle updates
   */
  for (a = 0, l0 = 18; a < 8; ++a, ++l0) {
    is = olim26::line3tris[a];
    if (nb[l0]) {
      for (b = 3, l1 = is[b]; b < 6; l1 = is[++b]) {
        if (nb[l1]) {
          T = min(T, this->tri23(VAL(l0), VAL(l1), s, h));
        }
      }
    }
  }

  /**
   * Degree (1, 2, 3) tetrahedron updates
   */
  for (a = 0, l0 = 18; a < 8; ++a, ++l0) {
    is = olim26::line3tris[a];
    if (nb[l0]) {
      for (b = 0, l1 = is[0], l2 = is[1];
           b < 6;
           ++b, l1 = is[b % 6], l2 = is[(b + 1) % 6]) {
        if (nb[l1] && nb[l2]) {
          if (b % 2 == 0) {
            T = min(T, this->tetra123(VAL(l1), VAL(l2), VAL(l0), s, h));
          } else {
            T = min(T, this->tetra123(VAL(l2), VAL(l1), VAL(l0), s, h));
          }
        }
      }
    }
  }

#ifdef PRINT_UPDATES
  printf("olim26_rect::update_impl: T <- %g\n", T);
#endif
}

#endif // __OLIM26_RECT_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
