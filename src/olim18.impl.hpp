#ifndef __OLIM18_IMPL_HPP__
#define __OLIM18_IMPL_HPP__

// #define PRINT_UPDATES 1

#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.macros.hpp"
#include "olim18.defs.hpp"

// neighbor order:
//
// N, E, U, S, W, D, DS, DW, DE, UE, UN, DN, SW, SE, NE, NW, UW, US
// 0, 1, 2, 3, 4, 5, 6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17
//
// the order of the degree 2 neighbors is critically
// important for the hash function DEG2NB to work

template <class node, class update_rules>
int olim18<node, update_rules>::di[] = {
  1, 0, 0, -1, 0, 0, -1, 0, 0, 0, 1, 1, -1, -1, 1, 1, 0, -1
};

template <class node, class update_rules>
int olim18<node, update_rules>::dj[] = {
  0, 1, 0, 0, -1, 0, 0, -1, 1, 1, 0, 0, -1, 1, 1, -1, -1, 0
};

template <class node, class update_rules>
int olim18<node, update_rules>::dk[] = {
  0, 0, 1, 0, 0, -1, -1, -1, -1, 1, 1, -1, 0, 0, 0, 0, 1, 1
};

template <class node, class update_rules>
void olim18<node, update_rules>::get_valid_neighbors(int i, int j, int k,
                                                     abstract_node ** nb) {
  int a, b, c;
  for (int l = 0; l < 18; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->is_valid(a, b, c)) {
      nb[l] = &this->operator()(a, b, c);
    }
  }
}

template <class node, class update_rules>
void olim18<node, update_rules>::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();
  int k = static_cast<node *>(n)->get_k();

  for (int l = 0; l < 18; ++l) {
    this->stage(i + di[l], j + dj[l], k + dk[l]);
  }

  int a, b, c;
  for (int l = 0; l < 18; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      this->update(a, b, c);
    }
  }
}

// Directions lying on the "equator" of the stencil:
static int eqdirs[4] = {
  olim18_defs::N, olim18_defs::E, olim18_defs::S, olim18_defs::W};

// Look-up table for triangles incident on 1-lines:
// TODO: change these to use olim18_defs::dir enum
static int line1tris[6][4] = {
  {10, 11, 14, 15}, // N
  {8,  9,  13, 14}, // E
  {9,  10, 16, 17}, // U
  {6,  12, 13, 17}, // S
  {7,  12, 15, 16}, // W
  {6,  7,  8,  11}, // D
};

// Look-up table for triangles incident on 2-lines:
// TODO: change these to use olim18_defs::dir enum
static int line2tris[12][6] = {
  {3, 5, 7,  8,  12, 13}, // DS
  {4, 5, 6,  11, 12, 15}, // DW
  {1, 5, 6,  11, 13, 14}, // DE
  {1, 2, 10, 17, 13, 14}, // UE
  {0, 2, 9,  16, 14, 15}, // UN
  {0, 5, 7,  8,  14, 15}, // DN
  {3, 4, 6,  7,  16, 17}, // SW
  {1, 3, 6,  8,  9,  17}, // SE
  {0, 1, 8,  9,  10, 11}, // NE
  {0, 4, 7,  10, 11, 16}, // NW
  {2, 4, 10, 17, 12, 15}, // UW
  {2, 3, 9,  16, 12, 13}, // US
};

/**
 * Hash a pair of directions in {N, E, U, S, W, D} uniquely into the
 * range [6, 18). This assumes that i < j. The result corresponds to
 * the index of the 2nd degree direction in the enum olim18_defs::dir.
 */
#define DEG2NB(i, j) ((5*i*i + 8*i*j + 4*i + 9*j) % 13 + 5)

template <class node, class update_rules>
void olim18<node, update_rules>::update_impl(int i, int j, int k, double & T) {
  using namespace olim18_defs;
  using std::min;
  using std::max;

#ifdef PRINT_UPDATES
  printf("olim18::update_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  abstract_node * nb[18];
  memset(nb, 0x0, 18*sizeof(abstract_node *));
  get_valid_neighbors(i, j, k, nb);

  double h = this->get_h(), s = this->speed(i, j, k), s_[18];
  int l, l0, l1, l2, lmin, lmax, l01, l02, l12;

  for (l = 0; l < 18; ++l) {
    if (nb[l]) {
      s_[l] = this->speed(i + di[l], j + dj[l], k + dk[l]);
    }
  }

  /**
   * line updates (degree 1 and 2)
   */
  for (l = 0; l < 6; ++l) {
    if (nb[l] && !(nb[line1tris[l][0]] || nb[line1tris[l][1]] ||
                   nb[line1tris[l][2]] || nb[line1tris[l][3]])) {
      T = min(T, this->line1(VAL(l), s, s_[l], h));
    }
  }
  for (l = 6, l0 = 0; l < 18; ++l, ++l0) {
    if (nb[l] && !(nb[line2tris[l0][0]] || nb[line2tris[l0][1]] ||
                   nb[line2tris[l0][2]] || nb[line2tris[l0][3]] ||
                   nb[line2tris[l0][4]] || nb[line2tris[l0][5]])) {
      T = min(T, this->line2(VAL(l), s, s_[l], h));
    }
  }
  
  /**
   * triangular updates ((1, 2) and (2, 2) 2-pt updates)
   *
   * TODO: find a simpler way to do this iteration
   */
  {
    int dirs[6] = {U, DEG2NB(N, U), N, DEG2NB(N, E), E, DEG2NB(E, U)};
    do_tri_updates(nb, dirs, s, s_, h, T);
    
    dirs[0] = U, dirs[1] = DEG2NB(U, W), dirs[2] = W, dirs[3] = DEG2NB(S, W),
      dirs[4] = S, dirs[5] = DEG2NB(U, S);
    do_tri_updates(nb, dirs, s, s_, h, T);

    dirs[0] = D, dirs[1] = DEG2NB(W, D), dirs[2] = W, dirs[3] = DEG2NB(N, W),
      dirs[4] = N, dirs[5] = DEG2NB(N, D);
    do_tri_updates(nb, dirs, s, s_, h, T);

    dirs[0] = D, dirs[1] = DEG2NB(S, D), dirs[2] = S, dirs[3] = DEG2NB(E, S),
      dirs[4] = E, dirs[5] = DEG2NB(E, D);
    do_tri_updates(nb, dirs, s, s_, h, T);
  }

  /**
   * tetrahedral updates
   */
  for (l0 = U; l0 != D; l0 = D) {
    for (l = 0, l1 = N, l2 = E;
         l < 4;
         ++l, l1 = eqdirs[l], l2 = eqdirs[(l + 1) % 4]) {
      lmin = min(l0, l1), lmax = max(l0, l1), l01 = DEG2NB(l0, l1);
      lmin = min(l0, l2), lmax = max(l0, l2), l02 = DEG2NB(l0, l2);
      lmin = min(l1, l2), lmax = max(l1, l2), l12 = DEG2NB(l1, l2);

      /*
       * (1, 2, 2) 3-pt updates
       */
      if (nb[l0] && nb[l01] && nb[l02]) {
        T = min(T, this->tetra122(
          VAL(l0), VAL(l01), VAL(l02), s, s_[l0], s_[l01], s_[l02], h));
      }
      if (nb[l12]) {
        if (nb[l1] && nb[l01]) {
          T = min(T, this->tetra122(
            VAL(l1), VAL(l01), VAL(l12), s, s_[l1], s_[l01], s_[l12], h));
        }
        if (nb[l2] && nb[l02]) {
          T = min(T, this->tetra122(
            VAL(l2), VAL(l02), VAL(l12), s, s_[l2], s_[l02], s_[l12], h));
        }
      }

      /*
       * (2, 2, 2) 3-pt update
       */
      if (nb[l01] && nb[l02] && nb[l12]) {
        T = min(T, this->tetra222(
          VAL(l01), VAL(l02), VAL(l12), s, s_[l01], s_[l02], s_[l12], h));
      }
    }
  }

#ifdef PRINT_UPDATES
  printf("olim18::update_impl: T <- %g\n", T);
#endif
}

/**
 * TODO: replace this with a macro once it's working
 */
template <class node, class update_rules>
void olim18<node, update_rules>::do_tri_updates(
  abstract_node const * const * nb, int const * dirs, double s,
  double const * s_, double h, double & T)
  const
{
  using std::min;
  int i, j, l0, l1;
  for (i = 0, j = 1; i < 6; ++i, j = (i + 1) % 6) {
    l0 = dirs[i], l1 = dirs[j];
    if (nb[l0] && nb[l1]) {
      if (i % 2 == 0) {
        T = min(T, this->tri12(VAL(l0), VAL(l1), s, s_[l0], s_[l1], h));
      } else {
        T = min(T, this->tri12(VAL(l1), VAL(l0), s, s_[l1], s_[l0], h));
      }
    }
  }
  for (i = 1, j = 3; i < 6; i += 2, j = (i + 2) % 6) {
    l0 = dirs[i], l1 = dirs[j];
    if (nb[l0] && nb[l1]) {
      T = min(T, this->tri22(VAL(l0), VAL(l1), s, s_[l0], s_[l1], h));
    }
  }
}

#endif // __OLIM18_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
