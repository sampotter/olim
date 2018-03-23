#include <gtest/gtest.h>

#include <src/config.hpp>

#include "common.defs.hpp"
#include "update_rules.tetra_updates.hpp"

static update_rules::rhr_tetra_updates updates;

#define P001 1
#define P010 2
#define P011 3
#define P100 4
#define P101 5
#define P110 6
#define P111 7

#define TETRA111(tetra_updates, u0, u1, u2, s, s0, s1, s2, h) \
  tetra_updates.tetra(                                        \
    u0, u1, u2, s, s0, s1, s2, h,                             \
    ffvec<P001> {}, ffvec<P010> {}, ffvec<P100> {})

#define TETRA122(tetra_updates, u0, u1, u2, s, s0, s1, s2, h) \
  tetra_updates.tetra(                                        \
    u0, u1, u2, s, s0, s1, s2, h,                             \
    ffvec<P001> {}, ffvec<P011> {}, ffvec<P101> {})

#define TETRA222(tetra_updates, u0, u1, u2, s, s0, s1, s2, h) \
  tetra_updates.tetra(                                        \
    u0, u1, u2, s, s0, s1, s2, h,                             \
    ffvec<P011> {}, ffvec<P101> {}, ffvec<P110> {})

#define TETRA123(tetra_updates, u0, u1, u2, s, s0, s1, s2, h) \
  tetra_updates.tetra(                                        \
    u0, u1, u2, s, s0, s1, s2, h,                             \
    ffvec<P001> {}, ffvec<P011> {}, ffvec<P111> {})

#if COLLECT_STATS
double get_T(update_return_t const & tmp) { return tmp.value; }
#else
double get_T(double tmp) { return tmp; }
#endif

template <class rules>
void tetra111_is_symmetric() {
  rules r;

  double u0 = 2.5453289254261224;
  double u1 = 2.5453289254261224;
  double u2 = 2.2844570503761732;
  double s = 1, s0 = 1, s1 = 1, s2 = 1;
  double h = 1;
  double val012 = get_T(TETRA111(r, u0, u1, u2, s, s0, s1, s2, h));
  double val120 = get_T(TETRA111(r, u1, u2, u0, s, s0, s1, s2, h));
  double val201 = get_T(TETRA111(r, u2, u0, u1, s, s0, s1, s2, h));
  double val210 = get_T(TETRA111(r, u2, u1, u0, s, s0, s1, s2, h));
  double val102 = get_T(TETRA111(r, u1, u0, u2, s, s0, s1, s2, h));
  double val021 = get_T(TETRA111(r, u0, u2, u1, s, s0, s1, s2, h));
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);

  u0 = 2.4;
  u1 = 2.2;
  u2 = 2.0;
  val012 = get_T(TETRA111(r, u0, u1, u2, s, s0, s1, s2, h));
  val120 = get_T(TETRA111(r, u1, u2, u0, s, s0, s1, s2, h));
  val201 = get_T(TETRA111(r, u2, u0, u1, s, s0, s1, s2, h));
  val210 = get_T(TETRA111(r, u2, u1, u0, s, s0, s1, s2, h));
  val102 = get_T(TETRA111(r, u1, u0, u2, s, s0, s1, s2, h));
  val021 = get_T(TETRA111(r, u0, u2, u1, s, s0, s1, s2, h));
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);
}

void tetra122_is_symmetric() {
  double u0 = 2.4;
  double u1 = 2.2;
  double u2 = 2.0;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;
  double val012 = get_T(TETRA122(updates, u0, u1, u2, s, s0, s1, s2, h));
  double val021 = get_T(TETRA122(updates, u0, u2, u1, s, s0, s1, s2, h));
  ASSERT_DOUBLE_EQ(val012, val021);
}

void tetra222_is_symmetric() {
  double u0 = 2.5453289254261224;
  double u1 = 2.5453289254261224;
  double u2 = 2.2844570503761732;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;
  double val012 = get_T(TETRA222(updates, u0, u1, u2, s, s0, s1, s2, h));
  double val120 = get_T(TETRA222(updates, u1, u2, u0, s, s0, s1, s2, h));
  double val201 = get_T(TETRA222(updates, u2, u0, u1, s, s0, s1, s2, h));
  double val210 = get_T(TETRA222(updates, u2, u1, u0, s, s0, s1, s2, h));
  double val102 = get_T(TETRA222(updates, u1, u0, u2, s, s0, s1, s2, h));
  double val021 = get_T(TETRA222(updates, u0, u2, u1, s, s0, s1, s2, h));
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);

  u0 = 2.4;
  u1 = 2.2;
  u2 = 2.0;
  val012 = get_T(TETRA222(updates, u0, u1, u2, s, s0, s1, s2, h));
  val120 = get_T(TETRA222(updates, u1, u2, u0, s, s0, s1, s2, h));
  val201 = get_T(TETRA222(updates, u2, u0, u1, s, s0, s1, s2, h));
  val210 = get_T(TETRA222(updates, u2, u1, u0, s, s0, s1, s2, h));
  val102 = get_T(TETRA222(updates, u1, u0, u2, s, s0, s1, s2, h));
  val021 = get_T(TETRA222(updates, u0, u2, u1, s, s0, s1, s2, h));
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);
}

void tetra111_works() {
  {
    double u0 = 1.867422661146497;
    double u1 = 1.872030476918322;
    double u2 = 1.874601455103048;
    double s = 1.0, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 1.0/7.0;
    double uhat = 1.95377665722661;
    ASSERT_DOUBLE_EQ(get_T(TETRA111(updates, u0, u1, u2, s, s0, s1, s2, h)), uhat);
  }
  {
    double u0 = 0.6714002359494359;
    double u1 = 0.667974148883546;
    double u2 = 0.652837863557498;
    double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 0.02040816326530612;
    double uhat = 0.6726606175825081;
    ASSERT_DOUBLE_EQ(get_T(TETRA111(updates, u0, u1, u2, s, s0, s1, s2, h)), uhat);
  }
  {
    double u0 = 0.8701508258299168;
    double u1 = 0.9034596879080383;
    double u2 = 0.8999244233472412;
    double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 0.03703703703703703;
    double uhat = 0.9064782064785435;
    ASSERT_DOUBLE_EQ(get_T(TETRA111(updates, u0, u1, u2, s, s0, s1, s2, h)), uhat);
  }
}

void tetra122_works() {
  update_rules::mp0_tetra_updates r;
  {
    double u0 = 1.9915638315627207;
    double u1 = 1.4142135623730949;
    double u2 = 2;
    double s = 1;
    double s0 = 1;
    double s1 = 1;
    double s2 = 1;
    double h = 1;
    double uhat = 2.781997898655415;
    ASSERT_DOUBLE_EQ(
      get_T(r.tetra(u0, u1, u2, s, s0, s1, s2, h,
                    ffvec<P001> {}, ffvec<P101> {}, ffvec<P110> {})),
      uhat);
  }
}

void tetra123_works() {
  ASSERT_DOUBLE_EQ(get_T(TETRA123(updates, 0, 0, 0, 1, 1, 1, 1, 1)), 1.0);
}

template <class rules>
void tetra222_works() {
  rules r;
  ASSERT_DOUBLE_EQ(get_T(TETRA222(r, 0, 0, 0, 1, 1, 1, 1, 1)), 2.0/sqrt(3));
  {
    double u0 = 0.18181818181818182;
    double u1 = 0.33244129540839806;
    double u2 = 0.33244129540839812;
    double h = 0.18181818181818182;
    double U = 0.4389479204314718;
    ASSERT_DOUBLE_EQ(get_T(TETRA222(r, u0, u1, u2, 1, 1, 1, 1, h)), U);
  }
}

update_rules::mp1_tetra_updates mp1;

void mp1_tetra123_works() {
  double u0, u1, u2, s, s0, s1, s2, h, U;

  u0 = 0.889163, u1 = 0.817579, u2 = 0.75, s = 1, s0 = 1, s1 = 1, s2 = 1,
    h = 0.25, U = 1.1189646747175703;
  ASSERT_DOUBLE_EQ(get_T(TETRA123(mp1, u0, u1, u2, s, s0, s1, s2, h)), U);
}

TEST (tetra_updates, tetra111_is_symmetric) {
  tetra111_is_symmetric<update_rules::mp0_tetra_updates>();
  tetra111_is_symmetric<update_rules::mp1_tetra_updates>();
  tetra111_is_symmetric<update_rules::rhr_tetra_updates>();
}

TEST (tetra_updates, tetra122_is_symmetric) {
  tetra122_is_symmetric();
}

TEST (tetra_updates, tetra222_is_symmetric) {
  tetra222_is_symmetric();
}

TEST (tetra_updates, tetra111_works) {
  tetra111_works();
}

TEST (tetra_updates, tetra122_works) {
  tetra122_works();
}

TEST (tetra_updates, tetra123_works) {
  tetra123_works();
}

TEST (tetra_updates, tetra222_works) {
  tetra222_works<update_rules::mp0_tetra_updates>();
  tetra222_works<update_rules::mp1_tetra_updates>();
  tetra222_works<update_rules::rhr_tetra_updates>();
}

TEST (tetra_updates, mp1_tetra123_works) {
  mp1_tetra123_works();
}
