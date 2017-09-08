#ifndef __OLIM18_IMPL_HPP__
#define __OLIM18_IMPL_HPP__

#include "common.macros.hpp"
#include "olim18.defs.hpp"

// neighbor order: U, UN, UE, US, UW, N, NE, E, SE, S, SW, W, NW, D,
// DN, DE, DS, DW

template <class node, class update_rules>
int olim18<node, update_rules>::di[] = {
  0, 1, 0, -1, 0, 1, 1, 0, -1, -1, -1, 0, 1, 0, 1, 0, -1, 0
};

template <class node, class update_rules>
int olim18<node, update_rules>::dj[] = {
  0, 0, 1, 0, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 0, 1, 0, -1
};

template <class node, class update_rules>
int olim18<node, update_rules>::dk[] = {
  1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1
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

template <class node, class update_rules>
void olim18<node, update_rules>::update_impl(int i, int j, int k, double & T) {
  using namespace olim18_defs;
  using std::min;

  abstract_node * nb[18];
  memset(nb, 0x0, 18*sizeof(abstract_node *));
  get_valid_neighbors(i, j, k, nb);

  double h = this->get_h(), s = this->speed(i, j, k);

  double s_[18];
  for (int l = 0; l < 18; ++l) {
    if (nb[l]) {
      s_[l] = this->speed(i + di[l], j + dj[l], k + dk[l]);
    }
  }

  bool nb12[24] = {
    nb[U] && nb[UN], nb[U] && nb[UE], nb[U] && nb[US], nb[U] && nb[UW],
    nb[N] && nb[UN], nb[N] && nb[NE], nb[N] && nb[NW], nb[N] && nb[DN],
    nb[E] && nb[UE], nb[E] && nb[NE], nb[E] && nb[SE], nb[E] && nb[DE],
    nb[S] && nb[US], nb[S] && nb[SE], nb[S] && nb[SW], nb[S] && nb[DS],
    nb[W] && nb[UW], nb[W] && nb[SW], nb[S] && nb[NW], nb[S] && nb[DW],
    nb[D] && nb[DN], nb[D] && nb[DE], nb[D] && nb[DS], nb[D] && nb[DW]
  };

  bool nb22[24] = {
    nb[UN] && nb[UE], nb[UE] && nb[US], nb[US] && nb[UW], nb[UW] && nb[UN],
    nb[UN] && nb[NE], nb[NE] && nb[DN], nb[DN] && nb[NW], nb[NW] && nb[UN],
    nb[UE] && nb[NE], nb[NE] && nb[DE], nb[DE] && nb[SE], nb[SE] && nb[UE],
    nb[US] && nb[SE], nb[SE] && nb[DS], nb[DS] && nb[SW], nb[SW] && nb[US],
    nb[UW] && nb[SW], nb[SW] && nb[DW], nb[DW] && nb[NW], nb[NW] && nb[UW],
    nb[DS] && nb[DE], nb[DE] && nb[DN], nb[DN] && nb[DW], nb[DW] && nb[DS]
  };

  bool nb122[24] = {
    nb[U] && nb[UN] && nb[UE], nb[U] && nb[UE] && nb[US],
    nb[U] && nb[US] && nb[UW], nb[U] && nb[UW] && nb[UW],
    nb[N] && nb[UN] && nb[NE], nb[N] && nb[NE] && nb[DN],
    nb[N] && nb[DN] && nb[NW], nb[N] && nb[NW] && nb[UN],
    nb[E] && nb[UE] && nb[NE], nb[E] && nb[NE] && nb[DE],
    nb[E] && nb[DE] && nb[SE], nb[E] && nb[SE] && nb[UE],
    nb[S] && nb[US] && nb[SE], nb[S] && nb[SE] && nb[DS],
    nb[S] && nb[DS] && nb[SW], nb[S] && nb[SW] && nb[US],
    nb[W] && nb[UW] && nb[SW], nb[W] && nb[SW] && nb[DW],
    nb[W] && nb[DW] && nb[NW], nb[W] && nb[NW] && nb[UW],
    nb[D] && nb[DS] && nb[DE], nb[D] && nb[DE] && nb[DN],
    nb[D] && nb[DN] && nb[DW], nb[D] && nb[DW] && nb[DS]
  };

  bool nb222[8] = {
    nb[UN] && nb[NE] && nb[UE], nb[UE] && nb[SE] && nb[US],
    nb[US] && nb[SW] && nb[UW], nb[UW] && nb[NW] && nb[UN],
    nb[DN] && nb[NW] && nb[DW], nb[DW] && nb[SW] && nb[DS],
    nb[DS] && nb[SE] && nb[DE], nb[DE] && nb[NE] && nb[DN]
  };

  // One point updates ('degree 1' nodes):

  if (!(nb12[U_UN] || nb12[U_UE] || nb12[U_US] || nb12[U_UW]) && nb[U]) {
    T = min(T, this->line1(VAL(U), s, s_[U], h));
  }
  if (!(nb12[N_UN] || nb12[N_NE] || nb12[N_NW] || nb12[N_DN]) && nb[N]) {
    T = min(T, this->line1(VAL(N), s, s_[N], h));
  }
  if (!(nb12[E_UE] || nb12[E_NE] || nb12[E_SE] || nb12[E_DE]) && nb[E]) {
    T = min(T, this->line1(VAL(E), s, s_[E], h));
  }
  if (!(nb12[S_US] || nb12[S_SE] || nb12[S_SW] || nb12[S_DS]) && nb[S]) {
    T = min(T, this->line1(VAL(S), s, s_[S], h));
  }
  if (!(nb12[W_UW] || nb12[W_SW] || nb12[W_NW] || nb12[W_DW]) && nb[W]) {
    T = min(T, this->line1(VAL(W), s, s_[W], h));
  }
  if (!(nb12[D_DN] || nb12[D_DE] || nb12[D_DS] || nb12[D_DW]) && nb[D]) {
    T = min(T, this->line1(VAL(D), s, s_[D], h));
  }

  // One point updates ('degree 2' nodes):

  if (!(nb12[U_UN] || nb12[N_UN]) && nb[UN]) {
    T = min(T, this->line2(VAL(UN), s, s_[UN], h));
  }
  if (!(nb12[U_UE] || nb12[E_UE]) && nb[UE]) {
    T = min(T, this->line2(VAL(UE), s, s_[UE], h));
  }
  if (!(nb12[U_US] || nb12[S_US]) && nb[US]) {
    T = min(T, this->line2(VAL(US), s, s_[US], h));
  }
  if (!(nb12[U_UW] || nb12[W_UW]) && nb[UW]) {
    T = min(T, this->line2(VAL(UW), s, s_[UW], h));
  }
  if (!(nb12[N_NE] || nb12[E_NE]) && nb[NE]) {
    T = min(T, this->line2(VAL(NE), s, s_[NE], h));
  }
  if (!(nb12[E_SE] || nb12[S_SE]) && nb[SE]) {
    T = min(T, this->line2(VAL(SE), s, s_[SE], h));
  }
  if (!(nb12[S_SW] || nb12[W_SW]) && nb[SW]) {
    T = min(T, this->line2(VAL(SW), s, s_[SW], h));
  }
  if (!(nb12[W_NW] || nb12[N_NW]) && nb[NW]) {
    T = min(T, this->line2(VAL(NW), s, s_[NW], h));
  }
  if (!(nb12[D_DN] || nb12[N_DN]) && nb[DN]) {
    T = min(T, this->line2(VAL(DN), s, s_[DN], h));
  }
  if (!(nb12[D_DE] || nb12[E_DE]) && nb[DE]) {
    T = min(T, this->line2(VAL(DE), s, s_[DE], h));
  }
  if (!(nb12[D_DS] || nb12[S_DS]) && nb[DS]) {
    T = min(T, this->line2(VAL(DS), s, s_[DS], h));
  }
  if (!(nb12[D_DW] || nb12[W_DW]) && nb[DW]) {
    T = min(T, this->line2(VAL(DW), s, s_[DW], h));
  }

  // Two point updates ('21' triangles):

  if (!nb122[U_UN_UE]) {
    if (!nb122[U_UW_UN] && nb12[U_UN]) {
      T = min(T, this->tri12(VAL(U), VAL(UN), s, s_[U], s_[UN], h));
    }
    if (!nb122[U_UE_US] && nb12[U_UE]) {
      T = min(T, this->tri12(VAL(U), VAL(UE), s, s_[U], s_[UE], h));
    }
  }
  if (!nb122[U_US_UW]) {
    if (!nb122[U_UE_US] && nb12[U_US]) {
      T = min(T, this->tri12(VAL(U), VAL(US), s, s_[U], s_[US], h));
    }
    if (!nb122[U_UW_UN] && nb12[U_UW]) {
      T = min(T, this->tri12(VAL(U), VAL(UW), s, s_[U], s_[UW], h));
    }
  }
  
  if (!nb122[N_UN_NE]) {
    if (!nb122[N_NW_UN] && nb12[N_UN]) {
      T = min(T, this->tri12(VAL(N), VAL(UN), s, s_[N], s_[UN], h));
    }
    if (!nb122[N_NE_DN] && nb12[N_NE]) {
      T = min(T, this->tri12(VAL(N), VAL(NE), s, s_[N], s_[NE], h));
    }
  }
  if (!nb122[N_DN_NW]) {
    if (!nb122[N_NE_DN] && nb12[N_DN]) {
      T = min(T, this->tri12(VAL(N), VAL(DN), s, s_[N], s_[DN], h));
    }
    if (!nb122[N_NW_UN] && nb12[N_NW]) {
      T = min(T, this->tri12(VAL(N), VAL(NW), s, s_[N], s_[NW], h));
    }
  }

  if (!nb122[E_UE_NE]) {
    if (!nb122[E_SE_UE] && nb12[E_UE]) {
      T = min(T, this->tri12(VAL(E), VAL(UE), s, s_[E], s_[UE], h));
    }
    if (!nb122[E_NE_DE] && nb12[E_NE]) {
      T = min(T, this->tri12(VAL(E), VAL(NE), s, s_[E], s_[NE], h));
    }
  }
  if (!nb122[E_DE_SE]) {
    if (!nb122[E_NE_DE] && nb12[E_DE]) {
      T = min(T, this->tri12(VAL(E), VAL(DE), s, s_[E], s_[DE], h));
    }
    if (!nb122[E_SE_UE] && nb12[E_SE]) {
      T = min(T, this->tri12(VAL(E), VAL(SE), s, s_[E], s_[SE], h));
    }
  }

  if (!nb122[S_US_SE]) {
    if (!nb122[S_SW_US] && nb12[S_US]) {
      T = min(T, this->tri12(VAL(S), VAL(US), s, s_[S], s_[US], h));
    }
    if (!nb122[S_SE_DS] && nb12[S_SE]) {
      T = min(T, this->tri12(VAL(S), VAL(SE), s, s_[S], s_[SE], h));
    }
  }
  if (!nb122[S_DS_SW]) {
    if (!nb122[S_SE_DS] && nb12[S_DS]) {
      T = min(T, this->tri12(VAL(S), VAL(DS), s, s_[S], s_[DS], h));
    }
    if (!nb122[S_SW_US] && nb12[S_SW]) {
      T = min(T, this->tri12(VAL(S), VAL(SW), s, s_[S], s_[SW], h));
    }
  }

  if (!nb122[W_UW_SW]) {
    if (!nb122[W_NW_UW] && nb12[W_UW]) {
      T = min(T, this->tri12(VAL(W), VAL(UW), s, s_[W], s_[UW], h));
    }
    if (!nb122[W_SW_DW] && nb12[W_SW]) {
      T = min(T, this->tri12(VAL(W), VAL(SW), s, s_[W], s_[SW], h));
    }
  }
  if (!nb122[W_DW_NW]) {
    if (!nb122[W_SW_DW] && nb12[W_DW]) {
      T = min(T, this->tri12(VAL(W), VAL(DW), s, s_[W], s_[DW], h));
    }
    if (!nb122[W_NW_UW] && nb12[W_NW]) {
      T = min(T, this->tri12(VAL(W), VAL(NW), s, s_[W], s_[NW], h));
    }
  }

  if (!nb122[D_DS_DE]) {
    if (!nb122[D_DW_DS] && nb12[D_DS]) {
      T = min(T, this->tri12(VAL(D), VAL(DS), s, s_[D], s_[DS], h));
    }
    if (!nb122[D_DE_DN] && nb12[D_DE]) {
      T = min(T, this->tri12(VAL(D), VAL(DE), s, s_[D], s_[DE], h));
    }
  }
  if (!nb122[D_DN_DW]) {
    if (!nb122[D_DE_DN] && nb12[D_DN]) {
      T = min(T, this->tri12(VAL(D), VAL(DN), s, s_[D], s_[DN], h));
    }
    if (!nb122[D_DW_DS] && nb12[D_DW]) {
      T = min(T, this->tri12(VAL(D), VAL(DW), s, s_[D], s_[DW], h));
    }
  }

  // Two point updates ('22' triangles):

  if (!nb222[UN_NE_UE]) {
    if (!nb[U] && nb22[UN_UE]) {
      T = min(T, this->tri22(VAL(UN), VAL(UE), s, s_[UN], s_[UE], h));
    }
    if (!nb[N] && nb22[UN_NE]) {
      T = min(T, this->tri22(VAL(UN), VAL(NE), s, s_[UN], s_[NE], h));
    }
    if (!nb[E] && nb22[UE_NE]) {
      T = min(T, this->tri22(VAL(NE), VAL(UE), s, s_[NE], s_[UE], h));
    }
  }

  if (!nb222[UE_SE_US]) {
    if (!nb[U] && nb22[UE_US]) {
      T = min(T, this->tri22(VAL(UE), VAL(US), s, s_[UE], s_[US], h));
    }
    if (!nb[E] && nb22[SE_UE]) {
      T = min(T, this->tri22(VAL(UE), VAL(SE), s, s_[UE], s_[SE], h));
    }
    if (!nb[S] && nb22[US_SE]) {
      T = min(T, this->tri22(VAL(SE), VAL(US), s, s_[SE], s_[US], h));
    }
  }

  if (!nb222[US_SW_UW]) {
    if (!nb[U] && nb22[US_UW]) {
      T = min(T, this->tri22(VAL(US), VAL(UW), s, s_[US], s_[UW], h));
    }
    if (!nb[S] && nb22[SW_US]) {
      T = min(T, this->tri22(VAL(US), VAL(SW), s, s_[US], s_[SW], h));
    }
    if (!nb[W] && nb22[UW_SW]) {
      T = min(T, this->tri22(VAL(SW), VAL(UW), s, s_[SW], s_[UW], h));
    }
  }
//
  if (!nb222[UW_NW_UN]) {
    if (!nb[U] && nb22[UW_UN]) {
      T = min(T, this->tri22(VAL(UW), VAL(UN), s, s_[UW], s_[UN], h));
    }
    if (!nb[W] && nb22[NW_UW]) {
      T = min(T, this->tri22(VAL(UW), VAL(NW), s, s_[UW], s_[NW], h));
    }
    if (!nb[N] && nb22[NW_UN]) {
      T = min(T, this->tri22(VAL(NW), VAL(UN), s, s_[NW], s_[UN], h));
    }
  }

  if (!nb222[DN_NW_DW]) {
    if (!nb[D] && nb22[DN_DW]) {
      T = min(T, this->tri22(VAL(DN), VAL(DW), s, s_[DN], s_[DW], h));
    }
    if (!nb[N] && nb22[DN_NW]) {
      T = min(T, this->tri22(VAL(DN), VAL(NW), s, s_[DN], s_[NW], h));
    }
    if (!nb[W] && nb22[DW_NW]) {
      T = min(T, this->tri22(VAL(NW), VAL(DW), s, s_[NW], s_[DW], h));
    }
  }

  if (!nb222[DW_SW_DS]) {
    if (!nb[D] && nb22[DW_DS]) {
      T = min(T, this->tri22(VAL(DW), VAL(DS), s, s_[DW], s_[DS], h));
    }
    if (!nb[W] && nb22[SW_DW]) {
      T = min(T, this->tri22(VAL(DW), VAL(SW), s, s_[DW], s_[SW], h));
    }
    if (!nb[S] && nb22[DS_SW]) {
      T = min(T, this->tri22(VAL(SW), VAL(DS), s, s_[SW], s_[DS], h));
    }
  }

  if (!nb222[DS_SE_DE]) {
    if (!nb[D] && nb22[DS_DE]) {
      T = min(T, this->tri22(VAL(DS), VAL(DE), s, s_[DS], s_[DE], h));
    }
    if (!nb[S] && nb22[SE_DS]) {
      T = min(T, this->tri22(VAL(DS), VAL(SE), s, s_[DS], s_[SE], h));
    }
    if (!nb[E] && nb22[DE_SE]) {
      T = min(T, this->tri22(VAL(SE), VAL(DE), s, s_[SE], s_[DE], h));
    }
  }

  if (!nb222[DE_NE_DN]) {
    if (!nb[D] && nb22[DE_DN]) {
      T = min(T, this->tri22(VAL(DE), VAL(DN), s, s_[DE], s_[DN], h));
    }
    if (!nb[E] && nb22[NE_DE]) {
      T = min(T, this->tri22(VAL(DE), VAL(NE), s, s_[DE], s_[NE], h));
    }
    if (!nb[N] && nb22[NE_DN]) {
      T = min(T, this->tri22(VAL(NE), VAL(DN), s, s_[DE], s_[DN], h));
    }
  }

  // Three point updates ('122' tetrahedra):

  if (nb122[U_UN_UE]) {
    T = min(T, this->tetra122(VAL(U), VAL(UN), VAL(UE), s, s_[U], s_[UN], s_[UE], h));
  }
  if (nb122[U_UE_US]) {
    T = min(T, this->tetra122(VAL(U), VAL(UE), VAL(US), s, s_[U], s_[UE], s_[US], h));
  }
  if (nb122[U_US_UW]) {
    T = min(T, this->tetra122(VAL(U), VAL(US), VAL(UW), s, s_[U], s_[US], s_[UW], h));
  }
  if (nb122[U_UW_UN]) {
    T = min(T, this->tetra122(VAL(U), VAL(UW), VAL(UN), s, s_[U], s_[UW], s_[UN], h));
  }
  if (nb122[N_UN_NE]) {
    T = min(T, this->tetra122(VAL(N), VAL(UN), VAL(NE), s, s_[N], s_[UN], s_[NE], h));
  }
  if (nb122[N_NE_DN]) {
    T = min(T, this->tetra122(VAL(N), VAL(NE), VAL(DN), s, s_[N], s_[NE], s_[DN], h));
  }
  if (nb122[N_DN_NW]) {
    T = min(T, this->tetra122(VAL(N), VAL(DN), VAL(NW), s, s_[N], s_[DN], s_[NW], h));
  }
  if (nb122[N_NW_UN]) {
    T = min(T, this->tetra122(VAL(N), VAL(NW), VAL(UN), s, s_[N], s_[NW], s_[UN], h));
  }
  if (nb122[E_UE_NE]) {
    T = min(T, this->tetra122(VAL(E), VAL(UE), VAL(NE), s, s_[E], s_[UE], s_[NE], h));
  }
  if (nb122[E_NE_DE]) {
    T = min(T, this->tetra122(VAL(E), VAL(NE), VAL(DE), s, s_[E], s_[NE], s_[DE], h));
  }
  if (nb122[E_DE_SE]) {
    T = min(T, this->tetra122(VAL(E), VAL(DE), VAL(SE), s, s_[E], s_[DE], s_[SE], h));
  }
  if (nb122[E_SE_UE]) {
    T = min(T, this->tetra122(VAL(E), VAL(SE), VAL(UE), s, s_[E], s_[SE], s_[UE], h));
  }
  if (nb122[S_US_SE]) {
    T = min(T, this->tetra122(VAL(S), VAL(US), VAL(SE), s, s_[S], s_[US], s_[SE], h));
  }
  if (nb122[S_SE_DS]) {
    T = min(T, this->tetra122(VAL(S), VAL(SE), VAL(DS), s, s_[S], s_[SE], s_[DS], h));
  }
  if (nb122[S_DS_SW]) {
    T = min(T, this->tetra122(VAL(S), VAL(DS), VAL(SW), s, s_[S], s_[DS], s_[SW], h));
  }
  if (nb122[S_SW_US]) {
    T = min(T, this->tetra122(VAL(S), VAL(SW), VAL(US), s, s_[S], s_[SW], s_[US], h));
  }
  if (nb122[W_UW_SW]) {
    T = min(T, this->tetra122(VAL(W), VAL(UW), VAL(SW), s, s_[W], s_[UW], s_[SW], h));
  }
  if (nb122[W_SW_DW]) {
    T = min(T, this->tetra122(VAL(W), VAL(SW), VAL(DW), s, s_[W], s_[SW], s_[DW], h));
  }
  if (nb122[W_DW_NW]) {
    T = min(T, this->tetra122(VAL(W), VAL(DW), VAL(NW), s, s_[W], s_[DW], s_[NW], h));
  }
  if (nb122[W_NW_UW]) {
    T = min(T, this->tetra122(VAL(W), VAL(NW), VAL(UW), s, s_[W], s_[NW], s_[UW], h));
  }
  if (nb122[D_DS_DE]) {
    T = min(T, this->tetra122(VAL(D), VAL(DS), VAL(DE), s, s_[D], s_[DS], s_[DE], h));
  }
  if (nb122[D_DE_DN]) {
    T = min(T, this->tetra122(VAL(D), VAL(DE), VAL(DN), s, s_[D], s_[DE], s_[DN], h));
  }
  if (nb122[D_DN_DW]) {
    T = min(T, this->tetra122(VAL(D), VAL(DN), VAL(DW), s, s_[D], s_[DN], s_[DW], h));
  }
  if (nb122[D_DW_DS]) {
    T = min(T, this->tetra122(VAL(D), VAL(DW), VAL(DS), s, s_[D], s_[DW], s_[DS], h));
  }

  // Three point updates ('222' tetrahedra):

  if (nb222[UN_NE_UE]) {
    T = min(T, this->tetra222(VAL(UN), VAL(NE), VAL(UE), s, s_[UN], s_[NE], s_[UE], h));
  }
  if (nb222[UE_SE_US]) {
    T = min(T, this->tetra222(VAL(UE), VAL(SE), VAL(US), s, s_[UE], s_[SE], s_[US], h));
  }
  if (nb222[US_SW_UW]) {
    T = min(T, this->tetra222(VAL(US), VAL(SW), VAL(UW), s, s_[US], s_[SW], s_[UW], h));
  }
  if (nb222[UW_NW_UN]) {
    T = min(T, this->tetra222(VAL(UW), VAL(NW), VAL(UN), s, s_[UW], s_[NW], s_[UN], h));
  }
  if (nb222[DN_NW_DW]) {
    T = min(T, this->tetra222(VAL(DN), VAL(NW), VAL(DW), s, s_[DN], s_[NW], s_[DW], h));
  }
  if (nb222[DW_SW_DS]) {
    T = min(T, this->tetra222(VAL(DW), VAL(SW), VAL(DS), s, s_[DW], s_[SW], s_[DS], h));
  }
  if (nb222[DS_SE_DE]) {
    T = min(T, this->tetra222(VAL(DS), VAL(SE), VAL(DE), s, s_[DS], s_[SE], s_[DE], h));
  }
  if (nb222[DE_NE_DN]) {
    T = min(T, this->tetra222(VAL(DE), VAL(NE), VAL(DN), s, s_[DE], s_[NE], s_[DN], h));
  }
}

#endif // __OLIM18_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
