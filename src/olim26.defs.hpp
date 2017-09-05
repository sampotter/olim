#ifndef __OLIM26_DEFS_HPP_HPP__
#define __OLIM26_DEFS_HPP_HPP__

namespace olim26_defs {
  // Keep these ordered by degree (= Hamming norm):
  enum dir {
    U, N, E, S, W, D, // degree 1
    UN, UE, US, UW, NE, SE, SW, NW, DN, DE, DS, DW, // degree 2
    UNE, USE, USW, UNW, DNE, DSE, DSW, DNW, // degree 3
  };

  constexpr int deg1start = 0, deg1end = 6, deg2start = deg1end,
    deg2end = 18, deg3start = deg2end, deg3end = 26;

  enum nb12 {
    N_NE, NE_E, E_SE, SE_S, S_SW, SW_W, W_NW, NW_N,
    U_UN, UN_N, N_DN, DN_D, D_DS, DS_S, S_US, US_U,
    U_UE, UE_E, E_DE, DE_D, D_DW, DW_W, W_UW, UW_U
  };

  enum nb13 {
    U_UNW, U_UNE, U_USE, U_USW,
    N_UNW, N_DNW, N_DNE, N_UNE,
    E_UNE, E_DNE, E_DSE, E_USE,
    S_USE, S_DSE, S_DSW, S_USW,
    W_USW, W_DSW, W_DNW, W_UNW,
    D_DSW, D_DSE, D_DNE, D_DNW
  };

  enum nb23 {
    UN_UNW, UN_UNE, UE_UNE, UE_USE,
    US_USE, US_USW, UW_USW, UW_UNW,
    NW_UNW, NE_UNE, SE_USE, SW_USW,
    NW_DNW, NE_DNE, SE_DSE, SW_DSW,
    DN_DNW, DN_DNE, DE_DNE, DE_DSE,
    DS_DSE, DS_DSW, DW_DSW, DW_DNW
  };

  enum nb123 {
    U_UN_UNE, U_UE_UNE, U_UE_USE, U_US_USE,
    U_US_USW, U_UW_USW, U_UW_UNW, U_UN_UNW,

    N_UN_UNE, E_UE_UNE, E_UE_USE, S_US_USE,
    S_US_USW, W_UW_USW, W_UW_UNW, N_UN_UNW,

    N_NE_UNE, E_NE_UNE, E_SE_USE, S_SE_USE,
    S_SW_USW, W_SW_USW, W_NW_UNW, N_NW_UNW,

    D_DN_DNE, D_DE_DNE, D_DE_DSE, D_DS_DSE,
    D_DS_DSW, D_DW_DSW, D_DW_DNW, D_DN_DNW,

    N_DN_DNE, E_DE_DNE, E_DE_DSE, S_DS_DSE,
    S_DS_DSW, W_DW_DSW, W_DW_DNW, N_DN_DNW,

    N_NE_DNE, E_NE_DNE, E_SE_DSE, S_SE_DSE,
    S_SW_DSW, W_SW_DSW, W_NW_DNW, N_NW_DNW
  };
}

#endif // __OLIM26_DEFS_HPP_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
