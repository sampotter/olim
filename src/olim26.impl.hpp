#ifndef __OLIM26_IMPL_HPP__
#define __OLIM26_IMPL_HPP__

#include <algorithm>

#include "common.macros.hpp"
#include "olim26.defs.hpp"

// neighbor order: U, UN, UNE, UE, USE, US, USW, UW, UNW, N, NE, E,
// SE, S, SW, W, NW, DN, DNE, DE, DSE, DS, DSW, DW, DNW, D

template <class node, class update_rules>
int olim26<node, update_rules>::di[] = {
  0, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1,
  -1, -1, 0, 1, 0
};

template <class node, class update_rules>
int olim26<node, update_rules>::dj[] = {
  0, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1,
  0, -1, -1, -1, 0
};

template <class node, class update_rules>
int olim26<node, update_rules>::dk[] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1,
  -1, -1, -1, -1, -1
};

template <class node, class update_rules>
void olim26<node, update_rules>::get_valid_neighbors(int i, int j, int k,
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
void olim26<node, update_rules>::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();
  int k = static_cast<node *>(n)->get_k();

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

template <class node, class update_rules>
void olim26<node, update_rules>::update_impl(int i, int j, int k, double & T) {
  using namespace olim26_defs;
  using std::min;

  abstract_node * nb[26];
  memset(nb, 0x0, 26*sizeof(abstract_node *));
  get_valid_neighbors(i, j, k, nb);

  double h = this->get_h(), s = this->speed(i, j, k);

  double s_[26];
  for (int l = 0; l < 26; ++l) {
    if (nb[l]) {
      s_[l] = this->speed(i + di[l], j + dj[l], k + dk[l]);
    }
  }

  bool nb12[24] = {
    nb[N_NE], nb[NE_E], nb[E_SE], nb[SE_S],
    nb[S_SW], nb[SW_W], nb[W_NW], nb[NW_N],
    nb[U_UN], nb[UN_N], nb[N_DN], nb[DN_D],
    nb[D_DS], nb[DS_S], nb[S_US], nb[US_U],
    nb[U_UE], nb[UE_E], nb[E_DE], nb[DE_D],
    nb[D_DW], nb[DW_W], nb[W_UW], nb[UW_U]
  };

  bool nb13[24] = {
    nb[U] && nb[UNW], nb[U] && nb[UNE], nb[U] && nb[USE], nb[U] && nb[USW],
    nb[N] && nb[UNW], nb[N] && nb[DNW], nb[N] && nb[DNE], nb[N] && nb[UNE],
    nb[E] && nb[UNE], nb[E] && nb[DNE], nb[E] && nb[DSE], nb[E] && nb[USE],
    nb[S] && nb[USE], nb[S] && nb[DSE], nb[S] && nb[DSW], nb[S] && nb[USW],
    nb[W] && nb[USW], nb[W] && nb[DSW], nb[W] && nb[DNW], nb[W] && nb[UNW],
    nb[D] && nb[DSW], nb[D] && nb[DSE], nb[D] && nb[DNE], nb[D] && nb[DNW]
  };

  bool nb23[24] = {
    nb[UN] && nb[UNW], nb[UN] && nb[UNE], nb[UE] && nb[UNE], nb[UE] && nb[USE],
    nb[US] && nb[USE], nb[US] && nb[USW], nb[UW] && nb[USW], nb[UW] && nb[UNW],
    nb[NW] && nb[UNW], nb[NE] && nb[UNE], nb[SE] && nb[USE], nb[SW] && nb[USW],
    nb[NW] && nb[DNW], nb[NE] && nb[DNE], nb[SE] && nb[DSE], nb[SW] && nb[DSW],
    nb[DN] && nb[DNW], nb[DN] && nb[DNE], nb[DE] && nb[DNE], nb[DE] && nb[DSE],
    nb[DS] && nb[DSE], nb[DS] && nb[DSW], nb[DW] && nb[DSW], nb[DW] && nb[DNW]
  };

  bool nb123[48] = {
    nb[U] && nb[UN] && nb[UNE], nb[U] && nb[UE] && nb[UNE],
    nb[U] && nb[UE] && nb[USE], nb[U] && nb[US] && nb[USE],
    nb[U] && nb[US] && nb[USW], nb[U] && nb[UW] && nb[USW],
    nb[U] && nb[UW] && nb[UNW], nb[U] && nb[UN] && nb[UNW],
    nb[N] && nb[UN] && nb[UNE], nb[E] && nb[UE] && nb[UNE],
    nb[E] && nb[UE] && nb[USE], nb[S] && nb[US] && nb[USE],
    nb[S] && nb[US] && nb[USW], nb[W] && nb[UW] && nb[USW],
    nb[W] && nb[UW] && nb[UNW], nb[N] && nb[UN] && nb[UNW],
    nb[N] && nb[NE] && nb[UNE], nb[E] && nb[NE] && nb[UNE],
    nb[E] && nb[SE] && nb[USE], nb[S] && nb[SE] && nb[USE],
    nb[S] && nb[SW] && nb[USW], nb[W] && nb[SW] && nb[USW],
    nb[W] && nb[NW] && nb[UNW], nb[N] && nb[NW] && nb[UNW],
    nb[D] && nb[DN] && nb[DNE], nb[D] && nb[DE] && nb[DNE],
    nb[D] && nb[DE] && nb[DSE], nb[D] && nb[DS] && nb[DSE],
    nb[D] && nb[DS] && nb[DSW], nb[D] && nb[DW] && nb[DSW],
    nb[D] && nb[DW] && nb[DNW], nb[D] && nb[DN] && nb[DNW],
    nb[N] && nb[DN] && nb[DNE], nb[E] && nb[DE] && nb[DNE],
    nb[E] && nb[DE] && nb[DSE], nb[S] && nb[DS] && nb[DSE],
    nb[S] && nb[DS] && nb[DSW], nb[W] && nb[DW] && nb[DSW],
    nb[W] && nb[DW] && nb[DNW], nb[N] && nb[DN] && nb[DNW],
    nb[N] && nb[NE] && nb[DNE], nb[E] && nb[NE] && nb[DNE],
    nb[E] && nb[SE] && nb[DSE], nb[S] && nb[SE] && nb[DSE],
    nb[S] && nb[SW] && nb[DSW], nb[W] && nb[SW] && nb[DSW],
    nb[W] && nb[NW] && nb[DNW], nb[N] && nb[NW] && nb[DNW]
  };

  // One point updates:

  // TODO: for now, just compute all the one-point updates that are
  // present---later, we'll try to do this like we do the rest of our
  // olim methods (skipping unnecessary updates), but this will
  // require quite a lot of checks, so we're deprioritizing this for
  // now

  for (int i = deg1start; i < deg1end; ++i) {
    if (nb[i]) {
      T = min(T, this->line1(VAL(i), s, s_[i], h));
    }
  }

  for (int i = deg2start; i < deg2end; ++i) {
    if (nb[i]) {
      T = min(T, this->line2(VAL(i), s, s_[i], h));
    }
  }

  for (int i = deg3start; i < deg3end; ++i) {
    if (nb[i]) {
      T = min(T, this->line3(VAL(i), s, s_[i], h));
    }
  }

  // Two point updates ('degree 12' triangles):

  if (!(nb123[N_NE_UNE] || nb123[N_NE_DNE]) && nb12[N_NE]) {
    T = min(T, this->tri12(VAL(N), VAL(NE), s, s_[N], s_[NE], h));
  }
  if (!(nb123[E_NE_UNE] || nb123[E_NE_DNE]) && nb12[NE_E]) {
    T = min(T, this->tri12(VAL(E), VAL(NE), s, s_[E], s_[NE], h));
  }
  if (!(nb123[E_SE_USE] || nb123[E_SE_DSE]) && nb12[E_SE]) {
    T = min(T, this->tri12(VAL(E), VAL(SE), s, s_[E], s_[SE], h));
  }
  if (!(nb123[S_SE_USE] || nb123[S_SE_DSE]) && nb12[SE_S]) {
    T = min(T, this->tri12(VAL(S), VAL(SE), s, s_[S], s_[SE], h));
  }
  if (!(nb123[S_SW_USW] || nb123[S_SW_DSW]) && nb12[S_SW]) {
    T = min(T, this->tri12(VAL(S), VAL(SW), s, s_[S], s_[SW], h));
  }
  if (!(nb123[W_SW_USW] || nb123[W_SW_DSW]) && nb12[SW_W]) {
    T = min(T, this->tri12(VAL(SW), VAL(W), s, s_[SW], s_[W], h));
  }
  if (!(nb123[W_NW_UNW] || nb123[W_NW_DNW]) && nb12[W_NW]) {
    T = min(T, this->tri12(VAL(W), VAL(NW), s, s_[W], s_[NW], h));
  }
  if (!(nb123[N_NW_UNW] || nb123[N_NW_DNW]) && nb12[NW_N]) {
    T = min(T, this->tri12(VAL(NW), VAL(N), s, s_[NW], s_[N], h));
  }
  if (!(nb123[N_UN_UNE] || nb123[U_UN_UNW]) && nb12[U_UN]) {
    T = min(T, this->tri12(VAL(U), VAL(UN), s, s_[U], s_[UN], h));
  }
  if (!(nb123[N_UN_UNE] || nb123[N_UN_UNW]) && nb12[UN_N]) {
    T = min(T, this->tri12(VAL(UN), VAL(N), s, s_[UN], s_[N], h));
  }
  if (!(nb123[N_DN_DNE] || nb123[N_DN_DNW]) && nb12[N_DN]) {
    T = min(T, this->tri12(VAL(N), VAL(DN), s, s_[N], s_[DN], h));
  }
  if (!(nb123[D_DN_DNE] || nb123[D_DN_DNW]) && nb12[DN_D]) {
    T = min(T, this->tri12(VAL(DN), VAL(D), s, s_[DN], s_[D], h));
  }
  if (!(nb123[D_DS_DSE] || nb123[D_DS_DSW]) && nb12[D_DS]) {
    T = min(T, this->tri12(VAL(D), VAL(DS), s, s_[D], s_[DS], h));
  }
  if (!(nb123[S_DS_DSE] || nb123[S_DS_DSW]) && nb12[DS_S]) {
    T = min(T, this->tri12(VAL(DS), VAL(S), s, s_[DS], s_[S], h));
  }
  if (!(nb123[S_US_USE] || nb123[S_US_USW]) && nb12[S_US]) {
    T = min(T, this->tri12(VAL(S), VAL(US), s, s_[S], s_[US], h));
  }
  if (!(nb123[U_US_USE] || nb123[U_US_USW]) && nb12[US_U]) {
    T = min(T, this->tri12(VAL(US), VAL(U), s, s_[US], s_[U], h));
  }
  if (!(nb123[U_UE_UNE] || nb123[U_UE_USE]) && nb12[U_UE]) {
    T = min(T, this->tri12(VAL(U), VAL(UE), s, s_[U], s_[UE], h));
  }
  if (!(nb123[E_UE_UNE] || nb123[E_UE_USE]) && nb12[UE_E]) {
    T = min(T, this->tri12(VAL(UE), VAL(E), s, s_[UE], s_[E], h));
  }
  if (!(nb123[E_DE_DNE] || nb123[E_DE_DSE]) && nb12[E_DE]) {
    T = min(T, this->tri12(VAL(E), VAL(DE), s, s_[E], s_[DE], h));
  }
  if (!(nb123[D_DE_DNE] || nb123[D_DE_DSE]) && nb12[DE_D]) {
    T = min(T, this->tri12(VAL(DE), VAL(D), s, s_[DE], s_[D], h));
  }
  if (!(nb123[D_DW_DNW] || nb123[D_DW_DSW]) && nb12[D_DW]) {
    T = min(T, this->tri12(VAL(D), VAL(DW), s, s_[D], s_[DW], h));
  }
  if (!(nb123[W_DW_DNW] || nb123[W_DW_DSW]) && nb12[DW_W]) {
    T = min(T, this->tri12(VAL(DW), VAL(W), s, s_[DW], s_[W], h));
  }
  if (!(nb123[W_UW_UNW] || nb123[W_UW_USW]) && nb12[W_UW]) {
    T = min(T, this->tri12(VAL(W), VAL(UW), s, s_[W], s_[UW], h));
  }
  if (!(nb123[U_UW_UNW] || nb123[U_UW_USW]) && nb12[UW_U]) {
    T = min(T, this->tri12(VAL(UW), VAL(U), s, s_[UW], s_[U], h));
  }

  // Two point updates ('degree 13' triangles):

  if (!(nb123[U_UW_UNW] || nb123[U_UN_UNW]) && nb13[U_UNW]) {
    T = min(T, this->tri13(VAL(U), VAL(UNW), s, s_[U], s_[UNW], h));
  }
  if (!(nb123[U_UN_UNE] || nb123[U_UE_UNE]) && nb13[U_UNE]) {
    T = min(T, this->tri13(VAL(U), VAL(UNE), s, s_[U], s_[UNE], h));
  }
  if (!(nb123[U_UE_USE] || nb123[U_US_USE]) && nb13[U_USE]) {
    T = min(T, this->tri13(VAL(U), VAL(USE), s, s_[U], s_[USE], h));
  }
  if (!(nb123[U_US_USW] || nb123[U_UW_USW]) && nb13[U_USW]) {
    T = min(T, this->tri13(VAL(U), VAL(USW), s, s_[U], s_[USW], h));
  }
  if (!(nb123[N_UN_UNW] || nb123[N_NW_UNW]) && nb13[N_UNW]) {
    T = min(T, this->tri13(VAL(N), VAL(UNW), s, s_[N], s_[UNW], h));
  }
  if (!(nb123[N_NW_DNW] || nb123[N_DN_DNW]) && nb13[N_DNW]) {
    T = min(T, this->tri13(VAL(N), VAL(DNW), s, s_[N], s_[DNW], h));
  }
  if (!(nb123[N_DN_DNE] || nb123[N_NE_DNE]) && nb13[N_DNE]) {
    T = min(T, this->tri13(VAL(N), VAL(DNE), s, s_[N], s_[DNE], h));
  }
  if (!(nb123[N_NE_UNE] || nb123[N_UN_UNE]) && nb13[N_UNE]) {
    T = min(T, this->tri13(VAL(N), VAL(UNE), s, s_[N], s_[UNE], h));
  }
  if (!(nb123[E_UE_UNE] || nb123[E_NE_UNE]) && nb13[E_UNE]) {
    T = min(T, this->tri13(VAL(E), VAL(UNE), s, s_[E], s_[UNE], h));
  }
  if (!(nb123[E_NE_DNE] || nb123[E_DE_DNE]) && nb13[E_DNE]) {
    T = min(T, this->tri13(VAL(E), VAL(DNE), s, s_[E], s_[DNE], h));
  }
  if (!(nb123[E_DE_DSE] || nb123[E_SE_DSE]) && nb13[E_DSE]) {
    T = min(T, this->tri13(VAL(E), VAL(DSE), s, s_[E], s_[DSE], h));
  }
  if (!(nb123[E_SE_USE] || nb123[E_UE_USE]) && nb13[E_USE]) {
    T = min(T, this->tri13(VAL(E), VAL(USE), s, s_[E], s_[USE], h));
  }
  if (!(nb123[S_US_USE] || nb123[S_SE_USE]) && nb13[S_USE]) {
    T = min(T, this->tri13(VAL(S), VAL(USE), s, s_[S], s_[USE], h));
  }
  if (!(nb123[S_SE_DSE] || nb123[S_DS_DSE]) && nb13[S_DSE]) {
    T = min(T, this->tri13(VAL(S), VAL(DSE), s, s_[S], s_[DSE], h));
  }
  if (!(nb123[S_DS_DSW] || nb123[S_SW_DSW]) && nb13[S_DSW]) {
    T = min(T, this->tri13(VAL(S), VAL(DSW), s, s_[S], s_[DSW], h));
  }
  if (!(nb123[S_SW_USW] || nb123[S_US_USW]) && nb13[S_USW]) {
    T = min(T, this->tri13(VAL(S), VAL(USW), s, s_[S], s_[USW], h));
  }
  if (!(nb123[W_UW_USW] || nb123[W_SW_USW]) && nb13[W_USW]) {
    T = min(T, this->tri13(VAL(W), VAL(USW), s, s_[W], s_[USW], h));
  }
  if (!(nb123[W_SW_DSW] || nb123[W_DW_DSW]) && nb13[W_DSW]) {
    T = min(T, this->tri13(VAL(W), VAL(DSW), s, s_[W], s_[DSW], h));
  }
  if (!(nb123[W_DW_DNW] || nb123[W_NW_DNW]) && nb13[W_DNW]) {
    T = min(T, this->tri13(VAL(W), VAL(DNW), s, s_[W], s_[DNW], h));
  }
  if (!(nb123[W_NW_UNW] || nb123[W_UW_UNW]) && nb13[W_UNW]) {
    T = min(T, this->tri13(VAL(W), VAL(UNW), s, s_[W], s_[UNW], h));
  }
  if (!(nb123[D_DW_DSW] || nb123[D_DS_DSW]) && nb13[D_DSW]) {
    T = min(T, this->tri13(VAL(D), VAL(DSW), s, s_[D], s_[DSW], h));
  }
  if (!(nb123[D_DS_DSE] || nb123[D_DE_DSE]) && nb13[D_DSE]) {
    T = min(T, this->tri13(VAL(D), VAL(DSE), s, s_[D], s_[DSE], h));
  }
  if (!(nb123[D_DE_DNE] || nb123[D_DN_DNE]) && nb13[D_DNE]) {
    T = min(T, this->tri13(VAL(D), VAL(DNE), s, s_[D], s_[DNE], h));
  }
  if (!(nb123[D_DN_DNW] || nb123[D_DW_DNW]) && nb13[D_DNW]) {
    T = min(T, this->tri13(VAL(D), VAL(DNW), s, s_[D], s_[DNW], h));
  }

  // Two point updates ('degree 23' triangles):

  if (!(nb123[U_UN_UNW] || nb123[N_UN_UNW]) && nb23[UN_UNW]) {
    T = min(T, this->tri23(VAL(UN), VAL(UNW), s, s_[UN], s_[UNW], h));
  }
  if (!(nb123[U_UN_UNE] || nb123[N_UN_UNE]) && nb23[UN_UNE]) {
    T = min(T, this->tri23(VAL(UN), VAL(UNE), s, s_[UN], s_[UNE], h));
  }
  if (!(nb123[U_UE_UNE] || nb123[E_UE_UNE]) && nb23[UE_UNE]) {
    T = min(T, this->tri23(VAL(UE), VAL(UNE), s, s_[UE], s_[UNE], h));
  }
  if (!(nb123[U_UE_USE] || nb123[E_UE_USE]) && nb23[UE_USE]) {
    T = min(T, this->tri23(VAL(UE), VAL(USE), s, s_[UE], s_[USE], h));
  }
  if (!(nb123[U_US_USE] || nb123[S_US_USE]) && nb23[US_USE]) {
    T = min(T, this->tri23(VAL(US), VAL(USE), s, s_[US], s_[USE], h));
  }
  if (!(nb123[U_US_USW] || nb123[S_US_USW]) && nb23[US_USW]) {
    T = min(T, this->tri23(VAL(US), VAL(USW), s, s_[US], s_[USW], h));
  }
  if (!(nb123[U_UW_USW] || nb123[W_UW_USW]) && nb23[UW_USW]) {
    T = min(T, this->tri23(VAL(UW), VAL(USW), s, s_[UW], s_[USW], h));
  }
  if (!(nb123[U_UW_UNW] || nb123[W_UW_UNW]) && nb23[UW_UNW]) {
    T = min(T, this->tri23(VAL(UW), VAL(UNW), s, s_[UW], s_[UNW], h));
  }
  if (!(nb123[N_NW_UNW] || nb123[W_NW_UNW]) && nb23[NW_UNW]) {
    T = min(T, this->tri23(VAL(NW), VAL(UNW), s, s_[NW], s_[UNW], h));
  }
  if (!(nb123[N_NE_UNE] || nb123[E_NE_UNE]) && nb23[NE_UNE]) {
    T = min(T, this->tri23(VAL(NE), VAL(UNE), s, s_[NE], s_[UNE], h));
  }
  if (!(nb123[S_SE_USE] || nb123[E_SE_USE]) && nb23[SE_USE]) {
    T = min(T, this->tri23(VAL(SE), VAL(USE), s, s_[SE], s_[USE], h));
  }
  if (!(nb123[S_SW_USW] || nb123[W_SW_USW]) && nb23[SW_USW]) {
    T = min(T, this->tri23(VAL(SW), VAL(USW), s, s_[SW], s_[USW], h));
  }
  if (!(nb123[N_NW_DNW] || nb123[W_NW_DNW]) && nb23[NW_DNW]) {
    T = min(T, this->tri23(VAL(NW), VAL(DNW), s, s_[NW], s_[DNW], h));
  }
  if (!(nb123[N_NE_DNE] || nb123[E_NE_DNE]) && nb23[NE_DNE]) {
    T = min(T, this->tri23(VAL(NE), VAL(DNE), s, s_[NE], s_[DNE], h));
  }
  if (!(nb123[S_SE_DSE] || nb123[E_SE_DSE]) && nb23[SE_DSE]) {
    T = min(T, this->tri23(VAL(SE), VAL(DSE), s, s_[SE], s_[DSE], h));
  }
  if (!(nb123[S_SW_DSW] || nb123[W_SW_DSW]) && nb23[SW_DSW]) {
    T = min(T, this->tri23(VAL(SW), VAL(DSW), s, s_[SW], s_[DSW], h));
  }
  if (!(nb123[D_DN_DNW] || nb123[N_DN_DNW]) && nb23[DN_DNW]) {
    T = min(T, this->tri23(VAL(DN), VAL(DNW), s, s_[DN], s_[DNW], h));
  }
  if (!(nb123[D_DN_DNE] || nb123[N_DN_DNE]) && nb23[DN_DNE]) {
    T = min(T, this->tri23(VAL(DN), VAL(DNE), s, s_[DN], s_[DNE], h));
  }
  if (!(nb123[D_DE_DNE] || nb123[E_DE_DNE]) && nb23[DE_DNE]) {
    T = min(T, this->tri23(VAL(DE), VAL(DNE), s, s_[DE], s_[DNE], h));
  }
  if (!(nb123[D_DE_DSE] || nb123[E_DE_DSE]) && nb23[DE_DSE]) {
    T = min(T, this->tri23(VAL(DE), VAL(DSE), s, s_[DE], s_[DSE], h));
  }
  if (!(nb123[D_DS_DSE] || nb123[S_DS_DSE]) && nb23[DS_DSE]) {
    T = min(T, this->tri23(VAL(DS), VAL(DSE), s, s_[DS], s_[DSE], h));
  }
  if (!(nb123[D_DS_DSW] || nb123[S_DS_DSW]) && nb23[DS_DSW]) {
    T = min(T, this->tri23(VAL(DS), VAL(DSW), s, s_[DS], s_[DSW], h));
  }
  if (!(nb123[D_DW_DSW] || nb123[W_DW_DSW]) && nb23[DW_DSW]) {
    T = min(T, this->tri23(VAL(DW), VAL(DSW), s, s_[DW], s_[DSW], h));
  }
  if (!(nb123[D_DW_DNW] || nb123[W_DW_DNW]) && nb23[DW_DNW]) {
    T = min(T, this->tri23(VAL(DW), VAL(DNW), s, s_[DW], s_[DNW], h));
  }
  
  // Three point updates ('degree 123' tetrahedra):

  if (nb123[U_UN_UNE]) {
    T = min(T, this->tetra123(VAL(U), VAL(UN), VAL(UNE), s, s_[U], s_[UN], s_[UNE], h));
  }
  if (nb123[U_UE_UNE]) {
    T = min(T, this->tetra123(VAL(U), VAL(UE), VAL(UNE), s, s_[U], s_[UE], s_[UNE], h));
  }
  if (nb123[U_UE_USE]) {
    T = min(T, this->tetra123(VAL(U), VAL(UE), VAL(USE), s, s_[U], s_[UE], s_[USE], h));
  }
  if (nb123[U_US_USE]) {
    T = min(T, this->tetra123(VAL(U), VAL(US), VAL(USE), s, s_[U], s_[US], s_[USE], h));
  }
  if (nb123[U_US_USW]) {
    T = min(T, this->tetra123(VAL(U), VAL(US), VAL(USW), s, s_[U], s_[US], s_[USW], h));
  }
  if (nb123[U_UW_USW]) {
    T = min(T, this->tetra123(VAL(U), VAL(UW), VAL(USW), s, s_[U], s_[UW], s_[USW], h));
  }
  if (nb123[U_UW_UNW]) {
    T = min(T, this->tetra123(VAL(U), VAL(UW), VAL(UNW), s, s_[U], s_[UW], s_[UNW], h));
  }
  if (nb123[U_UN_UNW]) {
    T = min(T, this->tetra123(VAL(U), VAL(UN), VAL(UNW), s, s_[U], s_[UN], s_[UNW], h));
  }
  if (nb123[N_UN_UNE]) {
    T = min(T, this->tetra123(VAL(N), VAL(UN), VAL(UNE), s, s_[N], s_[UN], s_[UNE], h));
  }
  if (nb123[E_UE_UNE]) {
    T = min(T, this->tetra123(VAL(E), VAL(UE), VAL(UNE), s, s_[E], s_[UE], s_[UNE], h));
  }
  if (nb123[E_UE_USE]) {
    T = min(T, this->tetra123(VAL(E), VAL(UE), VAL(USE), s, s_[E], s_[UE], s_[USE], h));
  }
  if (nb123[S_US_USE]) {
    T = min(T, this->tetra123(VAL(S), VAL(US), VAL(USE), s, s_[S], s_[US], s_[USE], h));
  }
  if (nb123[S_US_USW]) {
    T = min(T, this->tetra123(VAL(S), VAL(US), VAL(USW), s, s_[S], s_[US], s_[USW], h));
  }
  if (nb123[W_UW_USW]) {
    T = min(T, this->tetra123(VAL(W), VAL(UW), VAL(USW), s, s_[W], s_[UW], s_[USW], h));
  }
  if (nb123[W_UW_UNW]) {
    T = min(T, this->tetra123(VAL(W), VAL(UW), VAL(UNW), s, s_[W], s_[UW], s_[UNW], h));
  }
  if (nb123[N_UN_UNW]) {
    T = min(T, this->tetra123(VAL(N), VAL(UN), VAL(UNW), s, s_[N], s_[UN], s_[UNW], h));
  }
  if (nb123[N_NE_UNE]) {
    T = min(T, this->tetra123(VAL(N), VAL(NE), VAL(UNE), s, s_[N], s_[NE], s_[UNE], h));
  }
  if (nb123[E_NE_UNE]) {
    T = min(T, this->tetra123(VAL(E), VAL(NE), VAL(UNE), s, s_[E], s_[NE], s_[UNE], h));
  }
  if (nb123[E_SE_USE]) {
    T = min(T, this->tetra123(VAL(E), VAL(SE), VAL(USE), s, s_[E], s_[SE], s_[USE], h));
  }
  if (nb123[S_SE_USE]) {
    T = min(T, this->tetra123(VAL(S), VAL(SE), VAL(USE), s, s_[S], s_[SE], s_[USE], h));
  }
  if (nb123[S_SW_USW]) {
    T = min(T, this->tetra123(VAL(S), VAL(SW), VAL(USW), s, s_[S], s_[SW], s_[USW], h));
  }
  if (nb123[W_SW_USW]) {
    T = min(T, this->tetra123(VAL(W), VAL(SW), VAL(USW), s, s_[W], s_[SW], s_[USW], h));
  }
  if (nb123[W_NW_UNW]) {
    T = min(T, this->tetra123(VAL(W), VAL(NW), VAL(UNW), s, s_[W], s_[NW], s_[UNW], h));
  }
  if (nb123[N_NW_UNW]) {
    T = min(T, this->tetra123(VAL(N), VAL(NW), VAL(UNW), s, s_[N], s_[NW], s_[UNW], h));
  }
  if (nb123[D_DN_DNE]) {
    T = min(T, this->tetra123(VAL(D), VAL(DN), VAL(DNE), s, s_[D], s_[DN], s_[DNE], h));
  }
  if (nb123[D_DE_DNE]) {
    T = min(T, this->tetra123(VAL(D), VAL(DE), VAL(DNE), s, s_[D], s_[DE], s_[DNE], h));
  }
  if (nb123[D_DE_DSE]) {
    T = min(T, this->tetra123(VAL(D), VAL(DE), VAL(DSE), s, s_[D], s_[DE], s_[DSE], h));
  }
  if (nb123[D_DS_DSE]) {
    T = min(T, this->tetra123(VAL(D), VAL(DS), VAL(DSE), s, s_[D], s_[DS], s_[DSE], h));
  }
  if (nb123[D_DS_DSW]) {
    T = min(T, this->tetra123(VAL(D), VAL(DS), VAL(DSW), s, s_[D], s_[DS], s_[DSW], h));
  }
  if (nb123[D_DW_DSW]) {
    T = min(T, this->tetra123(VAL(D), VAL(DW), VAL(DSW), s, s_[D], s_[DW], s_[DSW], h));
  }
  if (nb123[D_DW_DNW]) {
    T = min(T, this->tetra123(VAL(D), VAL(DW), VAL(DNW), s, s_[D], s_[DW], s_[DNW], h));
  }
  if (nb123[D_DN_DNW]) {
    T = min(T, this->tetra123(VAL(D), VAL(DN), VAL(DNW), s, s_[D], s_[DN], s_[DNW], h));
  }
  if (nb123[N_DN_DNE]) {
    T = min(T, this->tetra123(VAL(N), VAL(DN), VAL(DNE), s, s_[N], s_[DN], s_[DNE], h));
  }
  if (nb123[E_DE_DNE]) {
    T = min(T, this->tetra123(VAL(E), VAL(DE), VAL(DNE), s, s_[E], s_[DE], s_[DNE], h));
  }
  if (nb123[E_DE_DSE]) {
    T = min(T, this->tetra123(VAL(E), VAL(DE), VAL(DSE), s, s_[E], s_[DE], s_[DSE], h));
  }
  if (nb123[S_DS_DSE]) {
    T = min(T, this->tetra123(VAL(S), VAL(DS), VAL(DSE), s, s_[S], s_[DS], s_[DSE], h));
  }
  if (nb123[S_DS_DSW]) {
    T = min(T, this->tetra123(VAL(S), VAL(DS), VAL(DSW), s, s_[S], s_[DS], s_[DSW], h));
  }
  if (nb123[W_DW_DSW]) {
    T = min(T, this->tetra123(VAL(W), VAL(DW), VAL(DSW), s, s_[W], s_[DW], s_[DSW], h));
  }
  if (nb123[W_DW_DNW]) {
    T = min(T, this->tetra123(VAL(W), VAL(DW), VAL(DNW), s, s_[W], s_[DW], s_[DNW], h));
  }
  if (nb123[N_DN_DNW]) {
    T = min(T, this->tetra123(VAL(N), VAL(DN), VAL(DNW), s, s_[N], s_[DN], s_[DNW], h));
  }
  if (nb123[N_NE_DNE]) {
    T = min(T, this->tetra123(VAL(N), VAL(NE), VAL(DNE), s, s_[N], s_[NE], s_[DNE], h));
  }
  if (nb123[E_NE_DNE]) {
    T = min(T, this->tetra123(VAL(E), VAL(NE), VAL(DNE), s, s_[E], s_[NE], s_[DNE], h));
  }
  if (nb123[E_SE_DSE]) {
    T = min(T, this->tetra123(VAL(E), VAL(SE), VAL(DSE), s, s_[E], s_[SE], s_[DSE], h));
  }
  if (nb123[S_SE_DSE]) {
    T = min(T, this->tetra123(VAL(S), VAL(SE), VAL(DSE), s, s_[S], s_[SE], s_[DSE], h));
  }
  if (nb123[S_SW_DSW]) {
    T = min(T, this->tetra123(VAL(S), VAL(SW), VAL(DSW), s, s_[S], s_[SW], s_[DSW], h));
  }
  if (nb123[W_SW_DSW]) {
    T = min(T, this->tetra123(VAL(W), VAL(SW), VAL(DSW), s, s_[W], s_[SW], s_[DSW], h));
  }
  if (nb123[W_NW_DNW]) {
    T = min(T, this->tetra123(VAL(W), VAL(NW), VAL(DNW), s, s_[W], s_[NW], s_[DNW], h));
  }
  if (nb123[N_NW_DNW]) {
    T = min(T, this->tetra123(VAL(N), VAL(NW), VAL(DNW), s, s_[N], s_[NW], s_[DNW], h));
  }
}

#endif // __OLIM26_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
