#ifndef __OLIM6_DEFS_HPP__
#define __OLIM6_DEFS_HPP__

enum class dir {
  U, UN, UE, US, UW, N, NE, E, SE, S, SW, W, NW, D, DN, DE, DS, DW
};

enum class nb12 {
  U_UN, U_UE, U_US, U_UW,
  N_UN, N_NE, N_NW, N_DN,
  E_UE, E_NE, E_SE, E_DE,
  S_US, S_SE, S_SW, S_DS,
  W_UW, W_SW, S_NW, S_DW,
  D_DN, D_DE, D_DS, D_DW
};

enum class nb22 {
  UN_UE, UE_US, US_UW, UW_UN,
  UN_NE, NE_DN, DN_NW, NW_UN,
  UE_NE, NE_DE, DE_SE, SE_UE,
  US_SE, SE_DS, DS_SW, SW_US,
  UW_SW, SW_DW, DW_NW, NW_UW,
  DS_DE, DE_DN, DN_DW, DW_DS
};

enum class nb122 {
  U_UN_UE, U_UE_US, U_US_UW, U_UW_UW,
  N_UN_NE, N_NE_DN, N_DN_NW, N_NW_UN,
  E_UE_NE, E_NE_DE, E_DE_SE, E_SE_UE,
  S_US_SE, S_SE_DS, S_DS_SW, S_SW_US,
  W_UW_SW, W_SW_DW, W_DW_NW, W_NW_UW,
  D_DS_DE, D_DE_DN, D_DN_DW, D_DW_DS
};

enum class nb222 {
  UN_NE_UE,
  UE_SE_US,
  US_SW_UW,
  UW_NW_UN,
  DN_NW_DW,
  DW_SW_DS,
  DS_SE_DE,
  DE_NE_DN
};

#endif // __OLIM6_DEFS_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
