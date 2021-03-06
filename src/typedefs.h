#ifndef __TYPEDEFS_H__
#define __TYPEDEFS_H__

enum marcher_type {
  BASIC,
  OLIM4_MP0,
  OLIM4_RHR,
  OLIM4_RHR_LUT,
  OLIM6_MP0,
  OLIM6_MP1,
  OLIM6_RHR,
  OLIM8_MP0,
  OLIM8_MP1,
  OLIM8_RHR,
  OLIM18_MP0,
  OLIM18_MP1,
  OLIM18_RHR,
  OLIM26_MP0,
  OLIM26_MP1,
  OLIM26_RHR,
  SOLIM4_MP0
};

typedef double(* speed_func)(double, double);
typedef double(* speed_func_3d)(double, double, double);

#endif // __TYPEDEFS_H__
