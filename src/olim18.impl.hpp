#ifndef __OLIM18_RECT_IMPL_HPP__
#define __OLIM18_RECT_IMPL_HPP__

#include <src/config.hpp>

#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.macros.hpp"
#include "olim.macros.hpp"
#include "olim18.defs.hpp"

// neighbor order:
//
// N, E, U, S, W, D, DS, DW, DE, UE, UN, DN, SW, SE, NE, NW, UW, US
// 0, 1, 2, 3, 4, 5, 6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17
//
// the order of the degree 2 neighbors is critically
// important for the hash function DEG2NB to work

template <class node, class line_updates, class tri_updates, class tetra_updates>
int olim18_rect<node, line_updates, tri_updates, tetra_updates>::di[] = {
// N, E, U,  S, W, D
   1, 0, 0, -1, 0, 0,
// DS, DW, DE, UE, UN, DN, SW, SE, NE, NW, UW, US
   -1, 0,  0,  0,  1,  1, -1, -1,  1,  1,  0,  -1
};

template <class node, class line_updates, class tri_updates, class tetra_updates>
int olim18_rect<node, line_updates, tri_updates, tetra_updates>::dj[] = {
// N, E, U, S, W,  D
   0, 1, 0, 0, -1, 0,
// DS, DW, DE, UE, UN, DN, SW, SE, NE, NW, UW, US
   0,  -1, 1,  1,  0,  0,  -1, 1,  1,  -1, -1, 0
};

template <class node, class line_updates, class tri_updates, class tetra_updates>
int olim18_rect<node, line_updates, tri_updates, tetra_updates>::dk[] = {
// N, E, U, S, W, D
   0, 0, 1, 0, 0, -1,
// DS, DW, DE, UE, UN, DN, SW, SE, NE, NW, UW, US
   -1, -1, -1, 1,  1,  -1, 0,  0,  0,  0,  1,  1
};

template <class node, class line_updates, class tri_updates, class tetra_updates>
void olim18_rect<node, line_updates, tri_updates, tetra_updates>::get_valid_neighbors(
  int i, int j, int k, abstract_node ** nb)
{
  int a, b, c;
  for (int l = 0; l < 18; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->is_valid(a, b, c)) {
      nb[l] = &this->operator()(a, b, c);
    }
  }
}

template <class node, class line_updates, class tri_updates, class tetra_updates>
void olim18_rect<node, line_updates, tri_updates, tetra_updates>::stage_neighbors_impl(
  abstract_node * n)
{
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();
  int k = static_cast<node *>(n)->get_k();

#if PRINT_UPDATES
  printf("olim18_rect::stage_neighbors_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

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
static int eqdirs[4] = {olim18::N, olim18::E, olim18::S, olim18::W};

/**
 * Hash a pair of directions in {N, E, U, S, W, D} uniquely into the
 * range [6, 18). This assumes that i < j. The result corresponds to
 * the index of the 2nd degree direction in the enum olim18::dir.
 */
#define DEG2NB(i, j) ((5*i*i + 8*i*j + 4*i + 9*j) % 13 + 5)

template <class node, class line_updates, class tri_updates, class tetra_updates>
void olim18_rect<node, line_updates, tri_updates, tetra_updates>::update_impl(
  int i, int j, int k, double & T)
{
  using namespace olim18;
  using std::min;
  using std::max;

#ifdef PRINT_UPDATES
  printf("olim18_rect::update_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  abstract_node * nb[18];
  memset(nb, 0x0, 18*sizeof(abstract_node *));
  get_valid_neighbors(i, j, k, nb);

  double h = this->get_h(), s = this->speed(i, j, k), s_[18];
  for (int l = 0; l < 18; ++l) {
    if (nb[l]) {
      s_[l] = this->speed(i + di[l], j + dj[l], k + dk[l]);
    }
  }
  
  int l, l0, l1, l2, l01, l02, l12;

  /**
   * line updates (degree 1 and 2)
   */
  for (l = 0; l < 6; ++l) if (nb[l]) LINE1(l);
  for (l = 6, l0 = 0; l < 18; ++l, ++l0) if (nb[l]) LINE2(l);
  
  /**
   * (1, 2) triangular updates
   */
  {
    int dirs[8] = {U, UN, N, DN, D, DS, S, US};
    do_tri12_updates(nb, dirs, s_, s, h, T);

    dirs[1] = UE, dirs[2] = E, dirs[3] = DE, dirs[4] = D, dirs[5] = DW,
      dirs[6] = W, dirs[7] = UW;
    do_tri12_updates(nb, dirs, s_, s, h, T);

    dirs[0] = N, dirs[1] = NE, dirs[3] = SE, dirs[4] = S, dirs[5] = SW,
      dirs[7] = NW;
    do_tri12_updates(nb, dirs, s_, s, h, T);
  }

  /**
   * (2, 2) triangular updates
   */
  {
    int dirs[4] = {UN, UE, US, UW};
    do_tri22_updates(nb, dirs, s_, s, h, T);

    dirs[1] = NE, dirs[2] = DN, dirs[3] = NW;
    do_tri22_updates(nb, dirs, s_, s, h, T);

    dirs[0] = UW, dirs[1] = SW, dirs[2] = DW;
    do_tri22_updates(nb, dirs, s_, s, h, T);

    dirs[0] = US, dirs[2] = DS, dirs[3] = SE;
    do_tri22_updates(nb, dirs, s_, s, h, T);

    dirs[0] = UE, dirs[1] = NE, dirs[2] = DE;
    do_tri22_updates(nb, dirs, s_, s, h, T);

    dirs[0] = DW, dirs[1] = DN, dirs[3] = DS;
    do_tri22_updates(nb, dirs, s_, s, h, T);
  }

  // TODO: eliminate the code duplication below

  /**
   * upper octant tetrahedral updates
   */
  l0 = U;
  for (l = 0, l1 = N, l2 = E;
       l < 4;
       ++l, l1 = eqdirs[l], l2 = eqdirs[(l + 1) % 4]) {
    /**
     * Sort l0, l1, and l2 indices and hash them to obtain indices to
     * degree 2 neighbors.
     */
    l01 = DEG2NB(min(l0, l1), max(l0, l1));
    l02 = DEG2NB(min(l0, l2), max(l0, l2));
    l12 = DEG2NB(min(l1, l2), max(l1, l2));

    /*
     * (1, 2, 2) 3-pt updates
     */
    if (nb[l0] && nb[l01] && nb[l02]) TETRA122(l0, l01, l02);
    if (nb[l12]) {
      if (nb[l1] && nb[l01]) TETRA122(l1, l01, l12);
      if (nb[l2] && nb[l02]) TETRA122(l2, l02, l12);
    }

    /*
     * (2, 2, 2) 3-pt update
     */
    if (nb[l01] && nb[l02] && nb[l12]) TETRA222(l01, l02, l12);
  }

  /**
   * lower octant tetrahedral updates
   */
  l0 = D;
  for (l = 0, l1 = N, l2 = E;
       l < 4;
       ++l, l1 = eqdirs[l], l2 = eqdirs[(l + 1) % 4]) {
    /**
     * Sort l0, l1, and l2 indices and hash them to obtain indices to
     * degree 2 neighbors.
     */
    l01 = DEG2NB(min(l0, l1), max(l0, l1));
    l02 = DEG2NB(min(l0, l2), max(l0, l2));
    l12 = DEG2NB(min(l1, l2), max(l1, l2));

    /*
     * (1, 2, 2) 3-pt updates
     */
    if (nb[l0] && nb[l01] && nb[l02]) TETRA122(l0, l01, l02);
    if (nb[l12]) {
      if (nb[l1] && nb[l01]) TETRA122(l1, l01, l12);
      if (nb[l2] && nb[l02]) TETRA122(l2, l02, l12);
    }

    /*
     * (2, 2, 2) 3-pt update
     */
    if (nb[l01] && nb[l02] && nb[l12]) TETRA222(l01, l02, l12);
  }

#ifdef PRINT_UPDATES
  printf("olim18_rect::update_impl: T <- %g\n", T);
#endif
}

template <class node, class line_updates, class tri_updates, class tetra_updates>
void olim18_rect<node, line_updates, tri_updates, tetra_updates>::do_tri12_updates(
  abstract_node const * const * nb, int const * dirs, 
  double const * s_, double s, double h, double & T)
  const
{
  using std::min;
  int l0, l1;
  for (int i = 0, j = 1; i < 8; j = (++i + 1) % 8) {
    if (nb[l0 = dirs[i]] && nb[l1 = dirs[j]]) {
      if (i % 2 == 0) TRI12(l0, l1);
      else TRI12(l1, l0);
    }
  }
}

template <class node, class line_updates, class tri_updates, class tetra_updates>
void olim18_rect<node, line_updates, tri_updates, tetra_updates>::do_tri22_updates(
  abstract_node const * const * nb, int const * dirs,
  double const * s_, double s, double h, double & T)
  const
{
  using std::min;
  int l0, l1;
  for (int i = 0, j = 1; i < 4; j = (++i + 1) % 4) {
    if (nb[l0 = dirs[i]] && nb[l1 = dirs[j]]) TRI22(l0, l1);
  }
}

#endif // __OLIM18_RECT_IMPL_HPP__
