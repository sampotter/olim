#include <gtest/gtest.h>

#include <src/config.hpp>

#include "common.hpp"
#include "updates.tetra.hpp"

using namespace updates;

#define P001 1
#define P010 2
#define P011 3
#define P100 4
#define P101 5
#define P110 6
#define P111 7

// TODO: these aren't completely general... The order (as they are
// passed as template arguments) of the bit vectors may make a
// difference...

template <cost_func F>
info<2> tetra111(double u0, double u1, double u2, double s,
                 double s0, double s1, double s2, double h)
{
  F_wkspc<F, 2> w;
  set_args<F>(w, u0, u1, u2, s, s0, s1, s2, h);
  cost_functor_bv<F, 3, P001, P010, P100> func {w};
  updates::info<2> info;
  tetra_bv<F, 3, P001, P010, P100>()(func, info);
  return info;
}

template <cost_func F>
info<2> tetra122(double u0, double u1, double u2, double s,
                 double s0, double s1, double s2, double h)
{
  F_wkspc<F, 2> w;
  set_args<F>(w, u0, u1, u2, s, s0, s1, s2, h);
  cost_functor_bv<F, 3, P001, P011, P101> func {w};
  updates::info<2> info;
  tetra_bv<F, 3, P001, P011, P101>()(func, info);
  return info;
}

template <cost_func F>
info<2> tetra123(double u0, double u1, double u2, double s,
                 double s0, double s1, double s2, double h)
{
  F_wkspc<F, 2> w;
  set_args<F>(w, u0, u1, u2, s, s0, s1, s2, h);
  cost_functor_bv<F, 3, P001, P011, P111> func {w};
  updates::info<2> info;
  tetra_bv<F, 3, P001, P011, P111>()(func, info);
  return info;
}

template <cost_func F>
info<2> tetra222(double u0, double u1, double u2, double s,
                        double s0, double s1, double s2, double h)
{
  F_wkspc<F, 2> w;
  set_args<F>(w, u0, u1, u2, s, s0, s1, s2, h);
  cost_functor_bv<F, 3, P011, P101, P110> func {w};
  updates::info<2> info;
  tetra_bv<F, 3, P011, P101, P110>()(func, info);
  return info;
}

template <cost_func F>
info<2> tetra223(double u0, double u1, double u2, double s,
                 double s0, double s1, double s2, double h)
{
  F_wkspc<F, 2> w;
  set_args<F>(w, u0, u1, u2, s, s0, s1, s2, h);
  cost_functor_bv<F, 3, P011, P101, P111> func {w};
  updates::info<2> info;
  tetra_bv<F, 3, P011, P101, P111>()(func, info);
  return info;
}

TEST (updates_tetra, tetra122_works) {
  updates::info<2> info;
  double u0, u1, u2, s, s0, s1, s2, h;

  u0 = 0.23;
  u1 = 0.19;
  u2 = 0.33;
  s = 1;
  s0 = 1.05;
  s1 = 1.03;
  s2 = 1.01;
  h = 0.25;

  info = tetra122<MP0>(u0, u1, u2, s, s0, s1, s2, h);
  ASSERT_NEAR(info.lambda[0], 0.174095688817618, eps<double>);
  ASSERT_NEAR(info.lambda[1], -0.435239222044045, eps<double>);
  ASSERT_TRUE(isinf(info.value));

  info = tetra122<MP1>(u0, u1, u2, s, s0, s1, s2, h);
  ASSERT_NEAR(info.lambda[0], 0.168611759124756, 4*eps<double>);
  ASSERT_NEAR(info.lambda[1], 0.0, eps<double>);
  ASSERT_NEAR(info.value, 0.482695106395732, 4*eps<double>);

  info = tetra122<RHR>(u0, u1, u2, s, s0, s1, s2, h);
  ASSERT_NEAR(info.lambda[0], 0.177296882871175, 5*eps<double>);
  ASSERT_NEAR(info.lambda[1], -0.443242207177936, 3*eps<double>);
  ASSERT_TRUE(isinf(info.value));
}

TEST (updates_tetra, tetra222_works) {
  updates::info<2> info;
  double u0, u1, u2, s, s0, s1, s2, h;

  u0 = 0.95;
  u1 = 1.0;
  u2 = 0.9;
  s = 1;
  s0 = 1.05;
  s1 = 1.03;
  s2 = 1.01;
  h = 0.25;

  info = tetra222<MP0>(u0, u1, u2, s, s0, s1, s2, h);
  ASSERT_NEAR(info.lambda[0], 0.096421835877633, 4*eps<double>);
  ASSERT_NEAR(info.lambda[1], 0.570244830789034, eps<double>);
  ASSERT_NEAR(info.value, 1.231399034587304, eps<double>);

  info = tetra222<MP1>(u0, u1, u2, s, s0, s1, s2, h);
  ASSERT_NEAR(info.lambda[0], 0.095098593078116, 4*eps<double>);
  ASSERT_NEAR(info.lambda[1], 0.585931976496175, 5*eps<double>);
  ASSERT_NEAR(info.value, 1.230640841044546, eps<double>);

  info = tetra222<RHR>(u0, u1, u2, s, s0, s1, s2, h);
  ASSERT_NEAR(info.lambda[0], 0.092561627161795, 3*eps<double>);
  ASSERT_NEAR(info.lambda[1], 0.574105039504872, eps<double>);
  ASSERT_NEAR(info.value, 1.226887462097269, eps<double>);
}

template <cost_func F>
void tetra111_is_symmetric_with_constant_slowness() {
  double uhat = 3.022472858860726;
  double u0 = 2.5453289254261224;
  double u1 = 2.5453289254261224;
  double u2 = 2.2844570503761732;
  double s = 1, s0 = 1, s1 = 1, s2 = 1;
  double h = 1;
  double val012 = tetra111<F>(u0, u1, u2, s, s0, s1, s2, h).value;
  double val120 = tetra111<F>(u1, u2, u0, s, s0, s1, s2, h).value;
  double val201 = tetra111<F>(u2, u0, u1, s, s0, s1, s2, h).value;
  double val210 = tetra111<F>(u2, u1, u0, s, s0, s1, s2, h).value;
  double val102 = tetra111<F>(u1, u0, u2, s, s0, s1, s2, h).value;
  double val021 = tetra111<F>(u0, u2, u1, s, s0, s1, s2, h).value;
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);
  ASSERT_DOUBLE_EQ(val012, uhat);
  ASSERT_DOUBLE_EQ(val120, uhat);
  ASSERT_DOUBLE_EQ(val201, uhat);
  ASSERT_DOUBLE_EQ(val210, uhat);
  ASSERT_DOUBLE_EQ(val102, uhat);
  ASSERT_DOUBLE_EQ(val021, uhat);

  u0 = 2.4;
  u1 = 2.2;
  u2 = 2.0;
  val012 = tetra111<F>(u0, u1, u2, s, s0, s1, s2, h).value;
  val120 = tetra111<F>(u1, u2, u0, s, s0, s1, s2, h).value;
  val201 = tetra111<F>(u2, u0, u1, s, s0, s1, s2, h).value;
  val210 = tetra111<F>(u2, u1, u0, s, s0, s1, s2, h).value;
  val102 = tetra111<F>(u1, u0, u2, s, s0, s1, s2, h).value;
  val021 = tetra111<F>(u0, u2, u1, s, s0, s1, s2, h).value;
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);
}

TEST (updates_tetra, tetra111_is_symmetric_with_constant_slowness) {
  tetra111_is_symmetric_with_constant_slowness<MP0>();
  tetra111_is_symmetric_with_constant_slowness<MP1>();
  tetra111_is_symmetric_with_constant_slowness<RHR>();
}

template <cost_func F>
void tetra122_is_symmetric_with_constant_slowness() {
  double u0 = 2.4;
  double u1 = 2.2;
  double u2 = 2.0;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;
  double val012 = tetra122<F>(u0, u1, u2, s, s0, s1, s2, h).value;
  double val021 = tetra122<F>(u0, u2, u1, s, s0, s1, s2, h).value;
  ASSERT_DOUBLE_EQ(val012, val021);
}

TEST (updates_tetra, tetra122_is_symmetric_with_constant_slowness) {
  tetra122_is_symmetric_with_constant_slowness<MP0>();
  tetra122_is_symmetric_with_constant_slowness<MP1>();
  tetra122_is_symmetric_with_constant_slowness<RHR>();
}

template <cost_func F>
void tetra222_is_symmetric_with_constant_slowness() {
  double u0 = 2.5453289254261224;
  double u1 = 2.5453289254261224;
  double u2 = 2.2844570503761732;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;
  double val012 = tetra222<F>(u0, u1, u2, s, s0, s1, s2, h).value;
  double val120 = tetra222<F>(u1, u2, u0, s, s0, s1, s2, h).value;
  double val201 = tetra222<F>(u2, u0, u1, s, s0, s1, s2, h).value;
  double val210 = tetra222<F>(u2, u1, u0, s, s0, s1, s2, h).value;
  double val102 = tetra222<F>(u1, u0, u2, s, s0, s1, s2, h).value;
  double val021 = tetra222<F>(u0, u2, u1, s, s0, s1, s2, h).value;
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);

  u0 = 2.4;
  u1 = 2.2;
  u2 = 2.0;
  val012 = tetra222<F>(u0, u1, u2, s, s0, s1, s2, h).value;
  val120 = tetra222<F>(u1, u2, u0, s, s0, s1, s2, h).value;
  val201 = tetra222<F>(u2, u0, u1, s, s0, s1, s2, h).value;
  val210 = tetra222<F>(u2, u1, u0, s, s0, s1, s2, h).value;
  val102 = tetra222<F>(u1, u0, u2, s, s0, s1, s2, h).value;
  val021 = tetra222<F>(u0, u2, u1, s, s0, s1, s2, h).value;
  ASSERT_DOUBLE_EQ(val012, val120);
  ASSERT_DOUBLE_EQ(val120, val201);
  ASSERT_DOUBLE_EQ(val201, val210);
  ASSERT_DOUBLE_EQ(val210, val102);
  ASSERT_DOUBLE_EQ(val102, val021);
}

TEST (updates_tetra, tetra222_is_symmetric_with_constant_slowness) {
  tetra222_is_symmetric_with_constant_slowness<MP0>();
  tetra222_is_symmetric_with_constant_slowness<MP1>();
  tetra222_is_symmetric_with_constant_slowness<RHR>();
}

template <cost_func F>
void tetra223_is_symmetric_with_constant_slowness() {
  double u0 = 2.4;
  double u1 = 2.2;
  double u2 = 2.0;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;
  double val012 = tetra223<F>(u0, u1, u2, s, s0, s1, s2, h).value;
  double val102 = tetra223<F>(u1, u0, u2, s, s0, s1, s2, h).value;
  ASSERT_DOUBLE_EQ(val012, val102);
}

TEST (updates_tetra, tetra223_is_symmetric_with_constant_slowness) {
  tetra223_is_symmetric_with_constant_slowness<MP0>();
  tetra223_is_symmetric_with_constant_slowness<MP1>();
  tetra223_is_symmetric_with_constant_slowness<RHR>();
}

template <cost_func F>
void tetra111_works_with_constant_slowness() {
  {
    double u0 = 1.867422661146497;
    double u1 = 1.872030476918322;
    double u2 = 1.874601455103048;
    double s = 1.0, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 1.0/7.0;
    double uhat = 1.95377665722661;
    ASSERT_DOUBLE_EQ(tetra111<F>(u0, u1, u2, s, s0, s1, s2, h).value, uhat);
  }
  {
    double u0 = 0.6714002359494359;
    double u1 = 0.667974148883546;
    double u2 = 0.652837863557498;
    double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 0.02040816326530612;
    double uhat = 0.6726606175825081;
    ASSERT_DOUBLE_EQ(tetra111<F>(u0, u1, u2, s, s0, s1, s2, h).value, uhat);
  }
  {
    double u0 = 0.8701508258299168;
    double u1 = 0.9034596879080383;
    double u2 = 0.8999244233472412;
    double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 0.03703703703703703;
    double uhat = 0.9064782064785435;
    ASSERT_DOUBLE_EQ(tetra111<F>(u0, u1, u2, s, s0, s1, s2, h).value, uhat);
  }
}

TEST (updates_tetra, tetra111_works_with_constant_slowness) {
  tetra111_works_with_constant_slowness<MP0>();
  tetra111_works_with_constant_slowness<MP1>();
  tetra111_works_with_constant_slowness<RHR>();
}

template <cost_func F>
void tetra122_works_with_constant_slowness() {
  double u0 = 1.9915638315627207;
  double u1 = 1.4142135623730949;
  double u2 = 2;
  double s = 1;
  double s0 = 1;
  double s1 = 1;
  double s2 = 1;
  double h = 1;
  auto info = tetra122<F>(u0, u1, u2, s, s0, s1, s2, h);
  if (F == MP1) {
    ASSERT_NEAR(info.lambda[0], 0.707106781186548, 1e1*eps<double>);
    ASSERT_NEAR(info.lambda[1], 0.0, 1e1*eps<double>);
    ASSERT_NEAR(info.value, 2.8080604124904465, eps<double>);
  } else {
    ASSERT_NEAR(info.lambda[0], 0.707144527237664, 1e1*eps<double>);
    ASSERT_NEAR(info.lambda[1], -0.010332705568234, 1e1*eps<double>);
    ASSERT_TRUE(isinf(info.value));
  }
}

TEST (updates_tetra, tetra122_works_with_constant_slowness) {
  tetra122_works_with_constant_slowness<MP0>();
  tetra122_works_with_constant_slowness<MP1>();
  tetra122_works_with_constant_slowness<RHR>();
}

template <cost_func F>
void tetra123_works_with_constant_slowness() {
  double u = tetra123<F>(0, 0, 0, 1, 1, 1, 1, 1).value, uhat = 1.0;
  ASSERT_DOUBLE_EQ(u, uhat);
}

TEST (updates_tetra, tetra123_works_with_constant_slowness) {
  tetra123_works_with_constant_slowness<MP0>();
  tetra123_works_with_constant_slowness<MP1>();
  tetra123_works_with_constant_slowness<RHR>();
}

template <cost_func F>
void tetra222_works_with_constant_slowness() {
  {
    auto info = tetra222<F>(0, 0, 0, 1, 1, 1, 1, 1);
    ASSERT_NEAR(info.lambda[0], 1./3, eps<double>);
    ASSERT_NEAR(info.lambda[1], 1./3, eps<double>);
    ASSERT_NEAR(info.value, 2.0/sqrt(3), eps<double>);
  }
  {
    double u0 = 0.19;
    double u1 = 0.2;
    double u2 = 0.21;
    double h = 1./3;
    auto info = tetra222<F>(u0, u1, u2, 1, 1, 1, 1, h);
    ASSERT_NEAR(info.lambda[0], 1./3, eps<double>);
    ASSERT_NEAR(info.lambda[1], 0.298661098115352, eps<double>);
    ASSERT_NEAR(info.value, 0.584553613273210, eps<double>);
  }
}

TEST (updates_tetra, tetra222_works_with_constant_slowness) {
  tetra222_works_with_constant_slowness<MP0>();
  tetra222_works_with_constant_slowness<MP1>();
  tetra222_works_with_constant_slowness<RHR>();
}

template <cost_func F>
void tetra223_works_with_constant_slowness() {
  double u0 = 1.0;
  double u1 = 1.1;
  double u2 = 1.2;
  double h = 1./3;
  auto info = tetra223<F>(u0, u1, u2, 1, 1, 1, 1, h);
  if (F == MP1) {
    ASSERT_NEAR(info.lambda[0], 0.312009803085040, 1e1*eps<double>);
    ASSERT_NEAR(info.lambda[1], 0.0, eps<double>);
    ASSERT_NEAR(info.value, 1.448956973452861, eps<double>);
  } else {
    ASSERT_NEAR(info.lambda[0], 1.404519917477945, 1e1*eps<double>);
    ASSERT_NEAR(info.lambda[1], -2.213559752433835, 5*eps<double>);
    ASSERT_TRUE(isinf(info.value));
  }
}

TEST (updates_tetra, tetra223_works_with_constant_slowness) {
  tetra223_works_with_constant_slowness<MP0>();
  tetra223_works_with_constant_slowness<MP1>();
  tetra223_works_with_constant_slowness<RHR>();
}

template <cost_func F>
void tetra111_mp0_is_symmetric_with_nonconstant_slowness() {
  {
    double U0 = 0.1, U1 = 0.2, U2 = 0.3;
    double s = 1, s0 = 1.3, s1 = 1.2, s2 = 1.1;
    double h = 0.5;
    double u = 0.506865877325366;
    double Uhat012 = tetra111<F>(U0, U1, U2, s, s0, s1, s2, h).value;
    double Uhat120 = tetra111<F>(U1, U2, U0, s, s1, s2, s0, h).value;
    double Uhat201 = tetra111<F>(U2, U0, U1, s, s2, s0, s1, h).value;
    double Uhat021 = tetra111<F>(U0, U2, U1, s, s0, s2, s1, h).value;
    double Uhat210 = tetra111<F>(U2, U1, U0, s, s2, s1, s0, h).value;
    double Uhat102 = tetra111<F>(U1, U0, U2, s, s1, s0, s2, h).value;
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
    double u = 0.558883511494396;
    double Uhat012 = tetra111<F>(U0, U1, U2, s, s0, s1, s2, h).value;
    double Uhat120 = tetra111<F>(U1, U2, U0, s, s1, s2, s0, h).value;
    double Uhat201 = tetra111<F>(U2, U0, U1, s, s2, s0, s1, h).value;
    double Uhat021 = tetra111<F>(U0, U2, U1, s, s0, s2, s1, h).value;
    double Uhat210 = tetra111<F>(U2, U1, U0, s, s2, s1, s0, h).value;
    double Uhat102 = tetra111<F>(U1, U0, U2, s, s1, s0, s2, h).value;
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
    double u = 0.558883511494395;
    double Uhat012 = tetra111<F>(U0, U1, U2, s, s0, s1, s2, h).value;
    double Uhat120 = tetra111<F>(U1, U2, U0, s, s1, s2, s0, h).value;
    double Uhat201 = tetra111<F>(U2, U0, U1, s, s2, s0, s1, h).value;
    double Uhat021 = tetra111<F>(U0, U2, U1, s, s0, s2, s1, h).value;
    double Uhat210 = tetra111<F>(U2, U1, U0, s, s2, s1, s0, h).value;
    double Uhat102 = tetra111<F>(U1, U0, U2, s, s1, s0, s2, h).value;
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

template <cost_func F>
updates::info<2> do_tetra(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h)
{
  F_wkspc<F, 2> w;
  set_args<F>(w, u0, u1, u2, s, s0, s1, s2, h);

  geom_wkspc<2> g;
  g.init<3>(p0, p1, p2);

  cost_functor<F, 3, 2> func {w, g};

  if (F != MP1) {
    func.qr = new qr_wkspc<3, 2>();
    func.qr->init(p0, p1, p2);
  }

  updates::info<2> info;
  updates::tetra<F, 3>()(func, info);

  if (F != MP1) {
    delete func.qr;
  }

  return info;
}

TEST (updates_tetra, tetra111_is_symmetric_with_nonconstant_slowness) {
  tetra111_mp0_is_symmetric_with_nonconstant_slowness<MP0>();
  // TODO: add mp1 and rhr
}

template <cost_func F>
void tetra122_is_symmetric_with_nonconstant_slowness() {
  {
    double U0, U1, U2, s, s0, s1, s2, h, Uhat012, Uhat021;
    U0 = 0.5793323043350438;
    U1 = 0.5613607127519399;
    U2 = 0.5864779604924129;
    s = 0.05655157497705388;
    s0 = 0.1194445225749439;
    s1 = 0.1479946250737618;
    s2 = 0.0679609140327736;
    h = 0.2;
    Uhat012 = tetra122<F>(U0, U1, U2, s, s0, s1, s2, h).value;
    Uhat021 = tetra122<F>(U0, U2, U1, s, s0, s2, s1, h).value;
    ASSERT_EQ(Uhat012, Uhat021);
  }
  {
    double p0[3], p1[3], p2[3], U0, U1, U2, s, s0, s1, s2, h, Uhat1, Uhat2;

    p0[0] = 0;
    p0[1] = 1;
    p0[2] = 1;
    p1[0] = 0;
    p1[1] = 0;
    p1[2] = 1;
    p2[0] = -1;
    p2[1] = 0;
    p2[2] = 1;
    U0 = 0.5879149021639897;
    U1 = 0.5823885480524517;
    U2 = 0.5637199168161153;
    s = 0.05655157497705388;
    s0 = 0.05655157497705388;
    s1 = 0.0876808174859417;
    s2 = 0.169616336715433;
    h = 0.2;
    Uhat1 = do_tetra<F>(p0, p1, p2, U0, U1, U2, s, s0, s1, s2, h).value;

    p0[0] = 0;
    p0[1] = 1;
    p0[2] = 1;
    p1[0] = 0;
    p1[1] = 0;
    p1[2] = 1;
    p2[0] = 1;
    p2[1] = 0;
    p2[2] = 1;
    U0 = 0.5879149021639898;
    U1 = 0.5823885480524518;
    U2 = 0.5637199168161156;
    s = 0.05655157497705388;
    s0 = 0.05655157497705388;
    s1 = 0.0876808174859417;
    s2 = 0.169616336715433;
    h = 0.2;
    Uhat2 = do_tetra<F>(p0, p1, p2, U0, U1, U2, s, s0, s1, s2, h).value;

    ASSERT_DOUBLE_EQ(Uhat1, Uhat2);
  }
}

TEST (updates_tetra, tetra122_is_symmetric_with_nonconstant_slowness) {
  tetra122_is_symmetric_with_nonconstant_slowness<MP0>();
  tetra122_is_symmetric_with_nonconstant_slowness<MP1>();
  tetra122_is_symmetric_with_nonconstant_slowness<RHR>();
}

TEST (updates_tetra, tetra123_works) {
  double p0[3], p1[3], p2[3], u0, u1, u2, s, s0, s1, s2, h;
  {
    p0[0] = 1;
    p0[1] = 0;
    p0[2] = 0;
    p1[0] = 1;
    p1[1] = -1;
    p1[2] = 0;
    p2[0] = 1;
    p2[1] = -1;
    p2[2] = 1;
    u0 = 0.55707818493084527;
    u1 = 0.547922866397345;
    u2 = 0.53876754786384473;
    s = 1;
    s0 = 1;
    s1 = 1;
    s2 = 1;
    h = 0.040000000000000001;
    auto info = do_tetra<MP0>(p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    ASSERT_NEAR(info.lambda[0], 0.0, eps<double>);
    ASSERT_NEAR(info.lambda[1], 0.241906380530053, eps<double>);
    ASSERT_NEAR(info.value, 0.594924720903115, eps<double>);
  }
}

template <cost_func F>
updates::info<2>
do_tetra(double const * p0, double const * p1, double const * p2,
         double u0, double u1, double u2, double s,
         double s0, double s1, double s2, double h,
         double const * p_fac, double s_fac)
{
  geom_fac_wkspc<2> g;
  g.init<3>(p0, p1, p2, p_fac);

  F_fac_wkspc<F, 2> w;
  set_args<F>(w, g, u0, u1, u2, s, s0, s1, s2, h, s_fac);

  cost_functor_fac<F, 3, 2> func {w, g};

  updates::info<2> info;
  updates::tetra<F, 3>()(func, info);

  return info;
}

template <cost_func F>
testing::AssertionResult basic_factoring_works(double tol = eps<double>) {
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

  auto info = do_tetra<F>(p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, p_fac, s_fac);
  {
    double gt = 1./3;
    if (fabs(info.lambda[0] - gt) > tol*gt + tol) {
      return testing::AssertionFailure()
        << "|" << info.lambda[0] << " - 1/3| > " << tol*gt + tol;
    }
  }
  {
    double gt = 1./3;
    if (fabs(info.lambda[1] - gt) > tol*gt + tol) {
      return testing::AssertionFailure()
        << "|" << info.lambda[1] << " - 1/3| > " << tol*gt + tol;
    }
  }
  {
    double gt = sqrt3;
    if (fabs(info.value - gt) > tol*gt + tol) {
      return testing::AssertionFailure()
        << "|" << info.value << " - sqrt(3)| > " << tol*gt + tol;
    }
  }
  return testing::AssertionSuccess();
}

TEST (updates_tetra, basic_factoring_test) {
  ASSERT_TRUE(basic_factoring_works<MP0>());
  ASSERT_TRUE(basic_factoring_works<MP1>());
  ASSERT_TRUE(basic_factoring_works<RHR>());
}
