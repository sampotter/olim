#include "common.defs.hpp"
#include "common.macros.hpp"
#include "test.hpp"
#include "update_rules.tri_updates.hpp"

update_rules::mp1_tri_updates mp1;
update_rules::rhr_tri_updates rhr;

#define P01 1
#define P10 2
#define P11 3

#define P001 1
#define P010 2
#define P011 3
#define P100 4
#define P101 5
#define P110 6
#define P111 7

#define TRI11(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P01> {}, ffvec<P10> {})

#define TRI12(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P01> {}, ffvec<P11> {})

#define TRI13(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P001> {}, ffvec<P111> {})

#define TRI22(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P011> {}, ffvec<P101> {})

#define TRI23(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P011> {}, ffvec<P111> {})

void rhr_tri11_is_symmetric() {
  double u0 = 0.0, u1 = 0.1, s = 1, s0 = 1, s1 = 1, h = 1;
  double t1 = TRI11(rhr, u0, u1, s, s0, s1, h);
  double t2 = TRI11(rhr, u1, u0, s, s1, s0, h);
  IS_APPROX_EQUAL(t1, t2);
}

void rhr_tri11_works() {
  IS_APPROX_EQUAL(TRI11(rhr, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.0);
  IS_APPROX_EQUAL(TRI11(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);
  IS_APPROX_EQUAL(TRI11(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.0);
}

void rhr_tri12_works() {
  double t;

  // TODO: add interior point test
  t = TRI12(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0);
  IS_APPROX_EQUAL(t, sqrt2);

  IS_APPROX_EQUAL(TRI12(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);

  double u0, u1, s, s0, s1, h;

  u0 = 0.923028, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(TRI12(rhr, u0, u1, s, s0, s1, h), 1.060660390593274);

  u0 = 0.852848, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(TRI12(rhr, u0, u1, s, s0, s1, h), 1.05597249);

  u0 = 0.923028, u1 = 0.817579, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(TRI12(rhr, u0, u1, s, s0, s1, h), 1.14970069);

  u0 = 0.883639, u1 = 0.817579, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(TRI12(rhr, u0, u1, s, s0, s1, h), 1.12475323);

  u0 = 0.852848, u1 = 0.65974, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(TRI12(rhr, u0, u1, s, s0, s1, h), 1.013293390593274);

  u0 = 0.883639, u1 = 0.65974, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(TRI12(rhr, u0, u1, s, s0, s1, h), 1.013293390593274);
}

void rhr_tri13_works() {
  // TODO: add interior point test
  IS_APPROX_EQUAL(TRI13(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.707106781186547);
  IS_APPROX_EQUAL(TRI13(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);
}

void rhr_tri22_works() {
  IS_APPROX_EQUAL(TRI22(rhr, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.414213562373095);
  IS_APPROX_EQUAL(TRI22(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.414213562373095);
  IS_APPROX_EQUAL(TRI22(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.414213562373095);

  double u0, u1, s, s0, s1, h;

  u0 = 0.65974, u1 = 0.817579, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(TRI22(rhr, u0, u1, s, s0, s1, h), 1.012639677310786);

  u0 = 0.65974, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(TRI22(rhr, u0, u1, s, s0, s1, h), 0.986849397845339);

  u0 = 0.707107, u1 = 0.817579, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(TRI22(rhr, u0, u1, s, s0, s1, h), 1.053198553345643);
}

void rhr_tri22_is_symmetric() {
  double u0 = 0.0, u1 = 0.1, s = 1, s0 = 1.0, s1 = 1.0, h = 1;
  double val01 = TRI22(rhr, u0, u1, s, s0, s1, h);
  double val10 = TRI22(rhr, u1, u0, s, s1, s0, h);
  IS_APPROX_EQUAL(val01, val10);
}

void rhr_tri23_works() {
  IS_APPROX_EQUAL(TRI23(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.732050807568877);
  IS_APPROX_EQUAL(TRI23(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.414213562373095);
}

void mp1_tri11_works() {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.1, u1 = 0, s = 1, s0 = 1.2, s1 = 1.1, h = 0.9, T = 0.7309362364283436;
  IS_APPROX_EQUAL(TRI11(mp1, u0, u1, s, s0, s1, h), T);
  IS_APPROX_EQUAL(TRI11(mp1, u0, u1, s, s0, s1, h), T);

  u0 = 0, u1 = 0.1, s = 1, s0 = 1.1, s1 = 1.3, h = 1.2, T = 0.97877243;
  IS_APPROX_EQUAL(TRI11(mp1, u0, u1, s, s0, s1, h), T);
  IS_APPROX_EQUAL(TRI11(mp1, u0, u1, s, s0, s1, h), T);

  u0 = 1.0, u1 = 0.853553, s = 1.0, s0 = 1.0, s1 = 1.0, h = 0.5, T = 1.27266;
  IS_APPROX_EQUAL(TRI11(mp1, u0, u1, s, s0, s1, h), T, 4e-6);
  IS_APPROX_EQUAL(TRI11(mp1, u0, u1, s, s0, s1, h), T, 4e-6);
}

void mp1_tri12_works() {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.46487865669951, u1 = 0.4, s = 1, s0 = 1, s1 = 1, h = 0.2;
  T = 0.6540631166978595;
  IS_APPROX_EQUAL(TRI12(mp1, u0, u1, s, s0, s1, h), T);
}

void mp1_tri13_works() {
  double u0, u1, s, s0, s1, h;

  u0 = 0, u1 = 1.41421, s = 1, s0 = 1, s1 = 1, h = 1;
  IS_APPROX_EQUAL(TRI13(mp1, u0, u1, s, s0, s1, h), 1.0);
}

void mp1_tri22_works() {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.65974, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  T = 0.986849397845339;
  IS_APPROX_EQUAL(TRI22(mp1, u0, u1, s, s0, s1, h), T);
}

int main() {
  rhr_tri11_is_symmetric();
  rhr_tri11_works();
  rhr_tri12_works();
  rhr_tri13_works();
  rhr_tri22_is_symmetric();
  rhr_tri22_works();
  rhr_tri23_works();
  mp1_tri11_works();
  mp1_tri12_works();
  mp1_tri13_works();
  mp1_tri22_works();
}
