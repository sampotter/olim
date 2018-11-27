#include <gtest/gtest.h>

#include "numopt.hpp"
#include "updates.common.hpp"

TEST (numopt, test_qpe_baryplex) {
  double G[3], c[2], x[2];
  {
    G[0] = 1.322471807186778;
    G[1] = 0.6280482242356766;
    G[2] = 1.035762733269118;
    c[0] = -0.9821315257790478;
    c[1] = 0.6125112981669493;
    qpe_bary<2, 0>(G, c, x);
    ASSERT_DOUBLE_EQ(x[0], 0.0);
    ASSERT_DOUBLE_EQ(x[1], -0.5913625567833622);
  }
  {
    G[0] = 1.47348599296532;
    G[1] = 0.2469229037436705;
    G[2] = 1.607389213768347;
    c[0] = -0.9930190065496006;
    c[1] = 0.9749502248113115;
    qpe_bary<2, 1>(G, c, x);
    ASSERT_DOUBLE_EQ(x[0], 0.67392497199868);
    ASSERT_DOUBLE_EQ(x[1], 0.0);
  }
  {
    G[0] = 1.915991244131425;
    G[1] = 0.2318001081857179;
    G[2] = 1.424349039815375;
    c[0] = -0.1131775740682571;
    c[1] = 1.635999657278292;
    qpe_bary<2, 2>(G, c, x);
    ASSERT_NEAR(x[0], 1.022590186764985, 1e-14);
    ASSERT_NEAR(x[1], -0.0225901867649847, 1e-14);
  }
}

TEST (numopt, test_qpi_baryplex) {
  bool error;
  double G[3], c[2], x[2], x0[2];
  for (int i = 0; i < 2; ++i) {
    if (i == 0) {
      x0[0] = x0[1] = 1./3;
    } else if (i == 1) {
      x0[0] = x0[1] = 0.0;
    }
    {
      G[0] = 0.2112628430553416;
      G[1] = -0.07940709555574356;
      G[2] = 0.3649291141863665;
      c[0] = -0.7577919164988186;
      c[1] = -0.5639689170274268;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_NEAR(x[0], 0.8682365591679995, 1e-14);
      ASSERT_NEAR(x[1], 0.1317634408320019, 1e-14);
    }
    {
      G[0] = 0.1766826028295286;
      G[1] = -0.007213173489582403;
      G[2] = 0.05466441540071013;
      c[0] = 0.3105083186004987;
      c[1] = -0.2494890644739389;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_NEAR(x[0], -4.440892098500626e-16, 1e-14);
      ASSERT_NEAR(x[1], 1.000000000000005, 1e-14);
    }
    {
      G[0] = 0.1398731795628713;
      G[1] = 0.1443894097074687;
      G[2] = 0.590026471661948;
      c[0] = -0.4335920223056836;
      c[1] = 0.3426244665386499;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_DOUBLE_EQ(x[0], 0.9999999999999999);
      ASSERT_DOUBLE_EQ(x[1], 0.0);
    }
    {
      G[0] = 0.7804349286530685;
      G[1] = 0.2283559117710858;
      G[2] = 0.6621076683127182;
      c[0] = 0.7147429038260958;
      c[1] = -0.2049660582997746;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_DOUBLE_EQ(x[0], 0.0);
      ASSERT_DOUBLE_EQ(x[1], 0.309566054146601);
    }
    {
      G[0] = 0.03998793788429245;
      G[1] = -0.05131392089981251;
      G[2] = 0.6514644398464839;
      c[0] = 0.7172386513288385;
      c[1] = 1.630235289164729;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_DOUBLE_EQ(x[0], 0.0);
      ASSERT_DOUBLE_EQ(x[1], 0.0);
    }
    {
      G[0] = 0.5996984878416955;
      G[1] = 0.07548854473125305;
      G[2] = 0.2674300462809748;
      c[0] = -0.1773751566188252;
      c[1] = -0.1960534878073328;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_DOUBLE_EQ(x[0], 0.2109897300222737);
      ASSERT_DOUBLE_EQ(x[1], 0.6735450359435763);
    }
    {
      G[0] = 0.3049276224231295;
      G[1] = -0.09291529749962908;
      G[2] = 0.7381116255079654;
      c[0] = 0.8350881650726819;
      c[1] = -0.2437151403779522;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_DOUBLE_EQ(x[0], 0.0);
      ASSERT_DOUBLE_EQ(x[1], 0.3301873753989831);
    }
    {
      G[0] = 0.9288870480671484;
      G[1] = 0.02769674842827917;
      G[2] = 0.7842908662640453;
      c[0] = -0.6668906707013855;
      c[1] = 0.1873310245789398;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_DOUBLE_EQ(x[0], 0.71794592473764);
      ASSERT_DOUBLE_EQ(x[1], 0.0);
    }
    {
      G[0] = 0.3120932443331088;
      G[1] = 0.02057818018813667;
      G[2] = 0.7934063383956034;
      c[0] = 0.1000928331393225;
      c[1] = -0.5445289299905477;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_DOUBLE_EQ(x[0], 0.0);
      ASSERT_DOUBLE_EQ(x[1], 0.6863178470336823);
    }
    {
      G[0] = 0.5111300221703463;
      G[1] = -0.1198364024041312;
      G[2] = 0.6875631691558609;
      c[0] = -2.138355269439939;
      c[1] = -0.8395887473366136;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_NEAR(x[0], 1.00000000000001, 1e-14);
      ASSERT_NEAR(x[1], -8.881784197001252e-16, 1e-14);
    }
    {
      G[0] = 0.3893554705283445;
      G[1] = 0.03674457678209533;
      G[2] = 0.2702975814212062;
      c[0] = 1.098424617888623;
      c[1] = -0.2778719327876389;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_NEAR(x[0], -4.440892098500626e-16, 1e-14);
      ASSERT_NEAR(x[1], 1.0, 1e-14);
    }
    {
      G[0] = 0.7925471686502944;
      G[1] = 0.2244581843371652;
      G[2] = 0.2128135976984587;
      c[0] = 0.2819840636705562;
      c[1] = 0.03347988224445142;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_DOUBLE_EQ(x[0], 0.0);
      ASSERT_DOUBLE_EQ(x[1], 0.0);
    }
    {
      G[0] = 0.4667997832499873;
      G[1] = -0.0549440710750627;
      G[2] = 0.4482577903690077;
      c[0] = -1.75021236844679;
      c[1] = -0.2856509715953298;
      qpi_bary<2>(G, c, x0, x, &error);
      ASSERT_FALSE(error);
      ASSERT_NEAR(x[0], 1.000000000000002, 1e-14);
      ASSERT_NEAR(x[1], 6.661338147750939e-16, 1e-14);
    }
    {
      G[0] = 1.8107622845933866;
      G[1] = 0.85973167304062859;
      G[2] = 1.7194633460812576;
    }
  }
}

TEST (numopt, sqp_bary_works_with_mp0) {
  bool error;
  updates::info<2> info;
  double u[3], h, s_hat, s[3], xgt[2], p[3][3], u_hat;
  {
    u[0] = 0.6463130101112646;
    u[1] = 0.962396190011243;
    u[2] = 0.9825909814628997;
    h = 0.4455862007108995;
    s_hat = 0.12299296536073;
    s[0] = 0.3028661333922572;
    s[1] = 0.2919026306839974;
    s[2] = 0.07245754527638229;
    xgt[0] = 0;
    xgt[1] = 0;
    u_hat = 0.7920957514401045;
    p[0][0] = -1.113500741486764;
    p[0][1] = -0.7696659137536819;
    p[0][2] = 1.117356138814467;
    p[1][0] = -0.006849328103348064;
    p[1][1] = 0.3713788127600577;
    p[1][2] = -1.089064295052236;
    p[2][0] = 1.53263030828475;
    p[2][1] = -0.2255844022712519;
    p[2][2] = 0.03255746416497347;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3>(w, p[0], p[1], p[2], u[0], u[1], u[2], s_hat, s[0], s[1], s[2], h);
    cost_functor<RHR, 3> func {w, p[0], p[1], p[2]};
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    assert(!error);
    func.eval(info.value);
    ASSERT_FALSE(error);
  }
}

TEST (numopt, sqp_bary_works_with_mp1) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_rhr) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_mp0_111) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_mp0_123) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_mp0_222) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_mp1_111) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_mp1_123) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_mp1_222) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_rhr_111) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_rhr_123) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_rhr_222) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_mp0_fac) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_mp1_fac) {
  ASSERT_TRUE(false); // TODO: implement me
}

TEST (numopt, sqp_bary_works_with_rhr_fac) {
  ASSERT_TRUE(false); // TODO: implement me
}
