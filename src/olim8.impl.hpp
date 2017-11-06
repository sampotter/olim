#ifndef __OLIM8_IMPL_HPP__
#define __OLIM8_IMPL_HPP__

#include <algorithm>
#include <cassert>

#include <src/config.hpp>

#include "common.macros.hpp"
#include "olim.macros.hpp"

template <class line_updates, class tri_updates>
void olim8<line_updates, tri_updates>::update_impl(int i, int j, double & T) {
  abstract_node * nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  abstract_node * x0 = 0x0;
  abstract_node * x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double s = speed(i, j), s0, s1, u0, u1, h = get_h();

  /*
   * First, do the adjacent and diagonal single point updates (and
   * only if the update triangles upon which the update edges are
   * incident are *not* present):
   */
  for (int k = 0; k < 8; k += 2) { // adjacent
    if ((x0 = nb[k]) && !nb[(k + 6) % 8] && !nb[(k + 7) % 8] &&
        !nb[(k + 1) % 8] && !nb[(k + 2) % 8]) {
      u0 = x0->get_value();
      s0 = speed(i + di[k], j + dj[k]);
      T = std::min(T, this->template line<1>(u0, s, s0, h));
    }
  }
  for (int k = 1; k < 8; k += 2) { // diagonal
    if ((x0 = nb[k]) && !nb[(k + 7) % 8] && !nb[(k + 1) % 8]) {
      u0 = x0->get_value();
      s0 = speed(i + di[k], j + dj[k]);
      T = std::min(T, this->template line<2>(u0, s, s0, h));
    }
  }

  /**
   * Next, do the diagonal triangle updates.
   */
  for (int k = 0, l = 1; k < 8; k += 2, l = k + 1) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s0 = speed(i + di[k], j + dj[k]);
      s1 = speed(i + di[l], j + dj[l]);
      T = std::min(T, this->tri12(u0, u1, s, s0, s1, h));
    }
  }
  for (int k = 1, l = 2; k < 8; k += 2, l = (k + 1) % 8) {
    if ((x0 = nb[l]) && (x1 = nb[k])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s0 = speed(i + di[l], j + dj[l]);
      s1 = speed(i + di[k], j + dj[k]);
      T = std::min(T, this->tri12(u0, u1, s, s0, s1, h));
    }
  }

#ifdef OLIM8_ADJ_UPDATES
  /**
   * Finally, do the adjacent triangle updates.
   */
  for (int k = 0, l = 2; k < 8; k += 2, l = (k + 2) % 8) {
    if ((x0 = nb[k]) && (x1 = nb[l])) {
      u0 = x0->get_value();
      u1 = x1->get_value();
      s0 = speed(i + di[k], j + dj[k]);
      s1 = speed(i + di[l], j + dj[l]);
      T = std::min(T, this->tri11(u0, u1, s, s0, s1, h));
    }
  }
#endif // OLIM8_ADJ_UPDATES
}

template <class line_updates, class tri_updates>
void olim8hu<line_updates, tri_updates>::update_impl(int i, int j, double & T) {
  abstract_node * nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  abstract_node * x0 = 0x0;
  abstract_node * x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double s = speed(i, j), h = get_h(), u0, s0, s1, Tnew;

  /**
   * Do every one-point update, and keep track of the index of the
   * minimizing point.
   */
  int argmin = -1;
  for (int k = 0; k < 8; ++k) {
    if ((x0 = nb[k])) {
      u0 = x0->get_value();
      s0 = speed(i + di[k], j + dj[k]);
      Tnew = k % 2 == 0 ? this->template line<1>(u0, s, s0, h) :
        this->template line <2>(u0, s, s0, h);
      if (Tnew < T) {
        T = Tnew;
        argmin = k;
      }
    }
  }
  assert(argmin != -1);

  /**
   * Do the two-point updates for the triangles which are adjacent to
   * the point that had the minimum one-point update.
   */
  int kprev = (argmin + 7) % 8, knext = (argmin + 1) % 8;
  u0 = nb[argmin]->get_value();
  s0 = speed(i + di[argmin], j + dj[argmin]);
  if (argmin % 2 == 0) {
    if ((x1 = nb[kprev])) {
      s1 = speed(i + di[kprev], j + dj[kprev]);
      T = std::min(T, this->tri12(u0, x1->get_value(), s, s0, s1, h));
    }
    if ((x1 = nb[knext])) {
      s1 = speed(i + di[knext], j + dj[knext]);
      T = std::min(T, this->tri12(u0, x1->get_value(), s, s0, s1, h));
    }
  } else {
    x1 = nb[kprev];
    abstract_node * x2 = nb[knext];
    if (x1) {
      s1 = speed(i + di[kprev], j + dj[kprev]);
      T = std::min(T, this->tri12(x1->get_value(), u0, s, s1, s0, h));
    }
    if (x2) {
      s1 = speed(i + di[knext], j + dj[knext]);
      T = std::min(T, this->tri12(x2->get_value(), u0, s, s1, s0, h));
    }
    if (x1 && x2) {
      s0 = speed(i + di[kprev], j + dj[kprev]);
      s1 = speed(i + di[knext], j + dj[knext]);
      T = std::min(T, this->tri12(x1->get_value(), x2->get_value(),
                                   s, s0, s1, h));
    }
  }
}

template <class line_updates, class tri_updates>
void olim8lut<line_updates, tri_updates>::update_impl(int i, int j, double & T) {
  using std::min;

  abstract_node * nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  get_valid_neighbors(i, j, nb);
  double s = speed(i, j), h = get_h(), s_[8];
  for (int k = 0; k < 8; ++k) {
    if (nb[k]) {
      s_[k] = this->speed(i + di[k], j + dj[k]);
    }
  }

  // implementing constrained algorithm for now
  for (int k = 0, l = 1, m = 2; k < 8; k += 2, l = k + 1, m = (l + 1) % 8) {
    switch (!!nb[k] + 2*!!nb[l] + 2*!!nb[m]) {
    case 2:
      LINE2(k);
      break;
    case 4:
      TRI12(l, m);
      break;
    case 3:
    case 6:
      LINE1(k);
      break;
    case 5: 
    case 7:
      TRI12(k, l);
      break;
    default:
      break;
    }
  }
}

#endif // __OLIM8_IMPL_HPP__
