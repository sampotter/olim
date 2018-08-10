#include <gtest/gtest.h>

#include <src/config.hpp>

#include "common.defs.hpp"
#include "update_rules.tetra_updates.hpp"

using namespace update_rules;

#define P001 1
#define P010 2
#define P011 3
#define P100 4
#define P101 5
#define P110 6
#define P111 7

// TODO: these aren't completely general... The order of the bit
// vectors may make a difference...

#define TETRA111(tetra_updates, u0, u1, u2, s, s0, s1, s2, h) \
  tetra_updates.tetra(                                        \
    u0, u1, u2, s, s0, s1, s2, h,                             \
    ffvec<P001> {}, ffvec<P010> {}, ffvec<P100> {})

#define TETRA122(tetra_updates, u0, u1, u2, s, s0, s1, s2, h) \
  tetra_updates.tetra(                                        \
    u0, u1, u2, s, s0, s1, s2, h,                             \
    ffvec<P001> {}, ffvec<P011> {}, ffvec<P101> {})

#define TETRA123(tetra_updates, u0, u1, u2, s, s0, s1, s2, h) \
  tetra_updates.tetra(                                        \
    u0, u1, u2, s, s0, s1, s2, h,                             \
    ffvec<P001> {}, ffvec<P011> {}, ffvec<P111> {})

#define TETRA222(tetra_updates, u0, u1, u2, s, s0, s1, s2, h) \
  tetra_updates.tetra(                                        \
    u0, u1, u2, s, s0, s1, s2, h,                             \
    ffvec<P011> {}, ffvec<P101> {}, ffvec<P110> {})

#define TETRA223(tetra_updates, u0, u1, u2, s, s0, s1, s2, h) \
  tetra_updates.tetra(                                        \
    u0, u1, u2, s, s0, s1, s2, h,                             \
    ffvec<P011> {}, ffvec<P101> {}, ffvec<P111> {})

template <int d>
double get_T(update_info<d> const & tmp) { return tmp.value; }

template <class rules>
void tetra111_is_symmetric_with_constant_slowness() {
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

TEST (tetra_updates, tetra111_is_symmetric_with_constant_slowness) {
  tetra111_is_symmetric_with_constant_slowness<mp0_tetra_updates>();
  tetra111_is_symmetric_with_constant_slowness<mp1_tetra_updates>();
  tetra111_is_symmetric_with_constant_slowness<rhr_tetra_updates>();
}

template <class rules>
void tetra122_is_symmetric_with_constant_slowness() {
  rules r;

  double u0 = 2.4;
  double u1 = 2.2;
  double u2 = 2.0;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;

  double val012 = get_T(TETRA122(r, u0, u1, u2, s, s0, s1, s2, h));
  double val021 = get_T(TETRA122(r, u0, u2, u1, s, s0, s1, s2, h));
  ASSERT_DOUBLE_EQ(val012, val021);
}

TEST (tetra_updates, tetra122_is_symmetric_with_constant_slowness) {
  tetra122_is_symmetric_with_constant_slowness<mp0_tetra_updates>();
  tetra122_is_symmetric_with_constant_slowness<mp1_tetra_updates>();
  tetra122_is_symmetric_with_constant_slowness<rhr_tetra_updates>();
}

template <class rules>
void tetra222_is_symmetric_with_constant_slowness() {
  rules r;

  double u0 = 2.5453289254261224;
  double u1 = 2.5453289254261224;
  double u2 = 2.2844570503761732;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;
  double val012 = get_T(TETRA222(r, u0, u1, u2, s, s0, s1, s2, h));
  double val120 = get_T(TETRA222(r, u1, u2, u0, s, s0, s1, s2, h));
  double val201 = get_T(TETRA222(r, u2, u0, u1, s, s0, s1, s2, h));
  double val210 = get_T(TETRA222(r, u2, u1, u0, s, s0, s1, s2, h));
  double val102 = get_T(TETRA222(r, u1, u0, u2, s, s0, s1, s2, h));
  double val021 = get_T(TETRA222(r, u0, u2, u1, s, s0, s1, s2, h));
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);

  u0 = 2.4;
  u1 = 2.2;
  u2 = 2.0;
  val012 = get_T(TETRA222(r, u0, u1, u2, s, s0, s1, s2, h));
  val120 = get_T(TETRA222(r, u1, u2, u0, s, s0, s1, s2, h));
  val201 = get_T(TETRA222(r, u2, u0, u1, s, s0, s1, s2, h));
  val210 = get_T(TETRA222(r, u2, u1, u0, s, s0, s1, s2, h));
  val102 = get_T(TETRA222(r, u1, u0, u2, s, s0, s1, s2, h));
  val021 = get_T(TETRA222(r, u0, u2, u1, s, s0, s1, s2, h));
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);
}

TEST (tetra_updates, tetra222_is_symmetric_with_constant_slowness) {
  tetra222_is_symmetric_with_constant_slowness<mp0_tetra_updates>();
  tetra222_is_symmetric_with_constant_slowness<mp1_tetra_updates>();
  tetra222_is_symmetric_with_constant_slowness<rhr_tetra_updates>();
}

template <class rules>
void tetra223_is_symmetric_with_constant_slowness() {
  rules r;

  double u0 = 2.4;
  double u1 = 2.2;
  double u2 = 2.0;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;
  double val012 = get_T(TETRA223(r, u0, u1, u2, s, s0, s1, s2, h));
  double val102 = get_T(TETRA223(r, u1, u0, u2, s, s0, s1, s2, h));
  ASSERT_DOUBLE_EQ(val012, val102);
}

TEST (tetra_updates, tetra223_is_symmetric_with_constant_slowness) {
  tetra223_is_symmetric_with_constant_slowness<mp0_tetra_updates>();
  tetra223_is_symmetric_with_constant_slowness<mp1_tetra_updates>();
  tetra223_is_symmetric_with_constant_slowness<rhr_tetra_updates>();
}

template <class rules>
void tetra111_works_with_constant_slowness() {
  rules r;
  {
    double u0 = 1.867422661146497;
    double u1 = 1.872030476918322;
    double u2 = 1.874601455103048;
    double s = 1.0, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 1.0/7.0;
    double uhat = 1.95377665722661;
    ASSERT_DOUBLE_EQ(get_T(TETRA111(r, u0, u1, u2, s, s0, s1, s2, h)), uhat);
  }
  {
    double u0 = 0.6714002359494359;
    double u1 = 0.667974148883546;
    double u2 = 0.652837863557498;
    double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 0.02040816326530612;
    double uhat = 0.6726606175825081;
    ASSERT_DOUBLE_EQ(get_T(TETRA111(r, u0, u1, u2, s, s0, s1, s2, h)), uhat);
  }
  {
    double u0 = 0.8701508258299168;
    double u1 = 0.9034596879080383;
    double u2 = 0.8999244233472412;
    double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 0.03703703703703703;
    double uhat = 0.9064782064785435;
    ASSERT_DOUBLE_EQ(get_T(TETRA111(r, u0, u1, u2, s, s0, s1, s2, h)), uhat);
  }
}

TEST (tetra_updates, tetra111_works_with_constant_slowness) {
  tetra111_works_with_constant_slowness<mp0_tetra_updates>();
  tetra111_works_with_constant_slowness<mp1_tetra_updates>();
  tetra111_works_with_constant_slowness<rhr_tetra_updates>();
}

template <class rules>
void tetra122_works_with_constant_slowness() {
  rules r;
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

TEST (tetra_updates, tetra122_works_with_constant_slowness) {
  tetra122_works_with_constant_slowness<mp0_tetra_updates>();
  tetra122_works_with_constant_slowness<mp1_tetra_updates>();
  tetra122_works_with_constant_slowness<rhr_tetra_updates>();
}

template <class rules>
void tetra123_works_with_constant_slowness() {
  rules r;
  ASSERT_DOUBLE_EQ(get_T(TETRA123(r, 0, 0, 0, 1, 1, 1, 1, 1)), 1.0);
}

TEST (tetra_updates, tetra123_works_with_constant_slowness) {
  tetra123_works_with_constant_slowness<mp0_tetra_updates>();
  tetra123_works_with_constant_slowness<mp1_tetra_updates>();
  tetra123_works_with_constant_slowness<rhr_tetra_updates>();
}

template <class rules>
void tetra222_works_with_constant_slowness() {
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

TEST (tetra_updates, tetra222_works_with_constant_slowness) {
  tetra222_works_with_constant_slowness<mp0_tetra_updates>();
  tetra222_works_with_constant_slowness<mp1_tetra_updates>();
  tetra222_works_with_constant_slowness<rhr_tetra_updates>();
}

template <class rules>
void tetra223_works_with_constant_slowness() {
  rules r;
  {
    double u0 = 0.18181818181818182;
    double u1 = 0.33244129540839806;
    double u2 = 0.33244129540839812;
    double h = 0.18181818181818182;
    double U = 0.4389479204314718;
    ASSERT_DOUBLE_EQ(get_T(TETRA223(r, u0, u1, u2, 1, 1, 1, 1, h)), U);
  }
}

TEST (tetra_updates, tetra223_works_with_constant_slowness) {
  tetra223_works_with_constant_slowness<mp0_tetra_updates>();
  tetra223_works_with_constant_slowness<mp1_tetra_updates>();
  tetra223_works_with_constant_slowness<rhr_tetra_updates>();
}

void tetra111_mp0_is_symmetric_with_nonconstant_slowness() {
  mp0_tetra_updates r;
  {
    double U0 = 0.1, U1 = 0.2, U2 = 0.3;
    double s = 1, s0 = 1.3, s1 = 1.2, s2 = 1.1;
    double h = 0.5;
    double u = 0.5101107150834783;
    double Uhat012 = get_T(TETRA111(r, U0, U1, U2, s, s0, s1, s2, h));
    double Uhat120 = get_T(TETRA111(r, U1, U2, U0, s, s1, s2, s0, h));
    double Uhat201 = get_T(TETRA111(r, U2, U0, U1, s, s2, s0, s1, h));
    double Uhat021 = get_T(TETRA111(r, U0, U2, U1, s, s0, s2, s1, h));
    double Uhat210 = get_T(TETRA111(r, U2, U1, U0, s, s2, s1, s0, h));
    double Uhat102 = get_T(TETRA111(r, U1, U0, U2, s, s1, s0, s2, h));
    ASSERT_DOUBLE_EQ(Uhat012, u);
    ASSERT_DOUBLE_EQ(Uhat120, u);
    ASSERT_DOUBLE_EQ(Uhat201, u);
    ASSERT_DOUBLE_EQ(Uhat021, u);
    ASSERT_DOUBLE_EQ(Uhat210, u);
    ASSERT_DOUBLE_EQ(Uhat102, u);
    ASSERT_DOUBLE_EQ(Uhat012, Uhat120);
    ASSERT_DOUBLE_EQ(Uhat120, Uhat201);
    ASSERT_DOUBLE_EQ(Uhat201, Uhat021);
    ASSERT_DOUBLE_EQ(Uhat021, Uhat210);
    ASSERT_DOUBLE_EQ(Uhat210, Uhat102);
  }
  {
    double U0 = 0.5514651482575854;
    double U1 = 0.5419072788623589;
    double U2 = 0.5415495962762169;
    double s = 0.2382400185837108;
    double s0 = 0.2420057484596495;
    double s1 = 0.2740574756206167;
    double s2 = 0.2914785686731952;
    double h = 0.1;
    double u = 0.5590785434269516;
    double Uhat012 = get_T(TETRA111(r, U0, U1, U2, s, s0, s1, s2, h));
    double Uhat120 = get_T(TETRA111(r, U1, U2, U0, s, s1, s2, s0, h));
    double Uhat201 = get_T(TETRA111(r, U2, U0, U1, s, s2, s0, s1, h));
    double Uhat021 = get_T(TETRA111(r, U0, U2, U1, s, s0, s2, s1, h));
    double Uhat210 = get_T(TETRA111(r, U2, U1, U0, s, s2, s1, s0, h));
    double Uhat102 = get_T(TETRA111(r, U1, U0, U2, s, s1, s0, s2, h));
    ASSERT_DOUBLE_EQ(Uhat012, u);
    ASSERT_DOUBLE_EQ(Uhat120, u);
    ASSERT_DOUBLE_EQ(Uhat201, u);
    ASSERT_DOUBLE_EQ(Uhat021, u);
    ASSERT_DOUBLE_EQ(Uhat210, u);
    ASSERT_DOUBLE_EQ(Uhat102, u);
    ASSERT_DOUBLE_EQ(Uhat012, Uhat120);
    ASSERT_DOUBLE_EQ(Uhat120, Uhat201);
    ASSERT_DOUBLE_EQ(Uhat201, Uhat021);
    ASSERT_DOUBLE_EQ(Uhat021, Uhat210);
    ASSERT_DOUBLE_EQ(Uhat210, Uhat102);
  }
  {
    double U0 = 0.5514651482575853;
    double U1 = 0.5419072788623587;
    double U2 = 0.5415495962762168;
    double s = 0.2382400185837107;
    double s0 = 0.2420057484596494;
    double s1 = 0.2740574756206167;
    double s2 = 0.2914785686731952;
    double h = 0.1;
    double u = 0.5590785434269516;
    double Uhat012 = get_T(TETRA111(r, U0, U1, U2, s, s0, s1, s2, h));
    double Uhat120 = get_T(TETRA111(r, U1, U2, U0, s, s1, s2, s0, h));
    double Uhat201 = get_T(TETRA111(r, U2, U0, U1, s, s2, s0, s1, h));
    double Uhat021 = get_T(TETRA111(r, U0, U2, U1, s, s0, s2, s1, h));
    double Uhat210 = get_T(TETRA111(r, U2, U1, U0, s, s2, s1, s0, h));
    double Uhat102 = get_T(TETRA111(r, U1, U0, U2, s, s1, s0, s2, h));
    ASSERT_DOUBLE_EQ(Uhat012, u);
    ASSERT_DOUBLE_EQ(Uhat120, u);
    ASSERT_DOUBLE_EQ(Uhat201, u);
    ASSERT_DOUBLE_EQ(Uhat021, u);
    ASSERT_DOUBLE_EQ(Uhat210, u);
    ASSERT_DOUBLE_EQ(Uhat102, u);
    ASSERT_DOUBLE_EQ(Uhat012, Uhat120);
    ASSERT_DOUBLE_EQ(Uhat120, Uhat201);
    ASSERT_DOUBLE_EQ(Uhat201, Uhat021);
    ASSERT_DOUBLE_EQ(Uhat021, Uhat210);
    ASSERT_DOUBLE_EQ(Uhat210, Uhat102);
  }
}

TEST (tetra_updates, tetra111_is_symmetric_with_nonconstant_slowness) {
  tetra111_mp0_is_symmetric_with_nonconstant_slowness();
}

template <class rules>
testing::AssertionResult basic_factoring_works(double tol = EPS(double)) {
  double u0, u1, u2, s, s0, s1, s2, h, p0[3], p1[3], p2[3], p_fac[3], s_fac;
  u0 = u1 = u2 = sqrt2;
  h = 1;
  s = s0 = s1 = s2 = s_fac = 1;

  // p_hat = (1, 1, 1), p0 = (1, 1, 0), p1 = (1, 0, 1),
  // p1 = (0, 1, 1), p_fac = (0, 0, 0)
  p0[0] = p1[1] = p2[2] = -1;
  p0[1] = p1[2] = p2[0] = 0;
  p0[2] = p1[0] = p2[1] = 0;
  p_fac[0] = p_fac[1] = p_fac[2] = -1;

  rules r;
  auto update = r.tetra(u0, u1, u2, s, s0, s1, s2, h, p0, p1, p2, p_fac, s_fac);
  {
    double gt = 1./3;
    if (fabs(update.lambda[0] - gt) > tol*gt + tol) {
      return testing::AssertionFailure()
        << "|" << update.lambda[0] << " - 1/3| > " << tol*gt + tol;
    }
  }
  {
    double gt = 1./3;
    if (fabs(update.lambda[1] - gt) > tol*gt + tol) {
      return testing::AssertionFailure()
        << "|" << update.lambda[1] << " - 1/3| > " << tol*gt + tol;
    }
  }
  {
    double gt = sqrt3;
    if (fabs(update.value - gt) > tol*gt + tol) {
      return testing::AssertionFailure()
        << "|" << update.value << " - sqrt(3)| > " << tol*gt + tol;
    }
  }
  return testing::AssertionSuccess();
}

TEST (tetra_updates, mp0_basic_factoring_test) {
  ASSERT_TRUE(basic_factoring_works<mp0_tetra_updates>());
}

TEST (tetra_updates, mp1_basic_factoring_test) {
  ASSERT_TRUE(basic_factoring_works<mp1_tetra_updates>());
}

TEST (tetra_updates, rhr_basic_factoring_test) {
  ASSERT_TRUE(basic_factoring_works<rhr_tetra_updates>());
}
