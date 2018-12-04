#include <gtest/gtest.h>

#include "common.defs.hpp"
#include "numopt.hpp"
#include "updates.common.hpp"

#define P001 1
#define P010 2
#define P011 3
#define P100 4
#define P101 5
#define P110 6
#define P111 7

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
  double u0, u1, u2, h, s, s0, s1, s2, p0[3], p1[3], p2[3];
  bool error;
  {
    u0 = 0.9508944153781352;
    u1 = 1.23310528676206;
    u2 = 0.9890460617169776;
    h = 0.6356613888613769;
    s = 0.8667498969993187;
    s0 = 0.6311887342690112;
    s1 = 0.8568953449704014;
    s2 = 1.26494521859783;
    p0[0] = 0.6524510729686149;
    p0[1] = 0.6049906419082594;
    p0[2] = 0.3872454314831349;
    p1[0] = 0.1421871592905041;
    p1[1] = 0.02513498571020312;
    p1[2] = 0.4211122537652413;
    p2[0] = 0.1841002894275112;
    p2[1] = 0.7257752674694531;
    p2[2] = 0.3703626865151981;
    F_wkspc<MP0, 2> w;
    set_args<MP0, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor<MP0, 3> func {w, p0, p1, p2};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0.5155525754868719, 2.22045e-16);
    ASSERT_NEAR(info.lambda[1], 0.4844474245131281, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_mp1) {
  double u0, u1, u2, h, s, s0, s1, s2, p0[3], p1[3], p2[3];
  bool error;
  {
    u0 = 0.5218856736612792;
    u1 = 0.5237900948295752;
    u2 = 0.5228817994319479;
    h = 0.005670468906829118;
    s = 0.2089466739931354;
    s0 = 0.905153559004464;
    s1 = 0.9089833436754959;
    s2 = 0.9078099933658579;
    p0[0] = 0.1040115747793785;
    p0[1] = 0.7455460737017173;
    p0[2] = 0.7362674555966385;
    p1[0] = 0.5618614252816374;
    p1[1] = 0.1841940975155265;
    p1[2] = 0.5972113503378549;
    p2[0] = 0.2999369900897889;
    p2[1] = 0.1341229328286824;
    p2[2] = 0.2126015333588429;
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor<MP1, 3> func {w, p0, p1, p2};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0, 2e-16);
    ASSERT_NEAR(info.lambda[1], 1, 2e-16);
  }
}

TEST (numopt, sqp_bary_works_with_rhr) {
  double u0, u1, u2, h, s, s0, s1, s2, p0[3], p1[3], p2[3];
  bool error;
  {
    u0 = 0.8944477555673932;
    u1 = 0.9612474890048042;
    u2 = 1.083854293536179;
    h = 0.4856516698980177;
    s = 0.9273562249981249;
    s0 = 0.9174938324161169;
    s1 = 1.264042242742724;
    s2 = 1.217790415332483;
    p0[0] = 0.9360273266897697;
    p0[1] = 0.1247740406604926;
    p0[2] = 0.7305853615057071;
    p1[0] = 0.6464774324258138;
    p1[1] = 0.833151985669295;
    p1[2] = 0.3982822282187755;
    p2[0] = 0.7498222093606359;
    p2[1] = 0.8352205104781305;
    p2[2] = 0.3224603973622594;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor<RHR, 3> func {w, p0, p1, p2};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0.3774084317899167, 3e-15);
    ASSERT_NEAR(info.lambda[1], 0, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_mp0_111) {
  double u0, u1, u2, h, s, s0, s1, s2;
  bool error;
  {
    u0 = 0.3480076597161135;
    u1 = 0.366528214414082;
    u2 = 0.4826058276823887;
    h = 0.1522340128629465;
    s = 0.09427839006414607;
    s0 = 0.9300406261074891;
    s1 = 0.990785037206068;
    s2 = 0.9372567408877572;
    F_wkspc<MP0, 2> w;
    set_args<MP0, 3, P100, P010, P001>(w, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor_bv<MP0, 3, P100, P010, P001> func {w};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0.4167015912580874, 4e-16);
    ASSERT_NEAR(info.lambda[1], 0, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_mp0_123) {
  double u0, u1, u2, h, s, s0, s1, s2;
  bool error;
  {
    u0 = 0.108016694136759;
    u1 = 0.577476739367107;
    u2 = 0.2380098353892041;
    h = 0.9080522031867692;
    s = 0.5593705724030039;
    s0 = 0.004579623947323475;
    s1 = 0.7007669019392007;
    s2 = 0.7752519068977054;
    F_wkspc<MP0, 2> w;
    set_args<MP0, 3, P100, P110, P111>(w, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor_bv<MP0, 3, P100, P110, P111> func {w};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], -2.886579864025407e-15, 3e-15);
    ASSERT_NEAR(info.lambda[1], -2.220446049250313e-16, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_mp0_222) {
  double u0, u1, u2, h, s, s0, s1, s2;
  bool error;
  {
    u0 = 0.3001844012139;
    u1 = 0.4783382851963944;
    u2 = 0.6700693493325152;
    h = 0.4438458367269573;
    s = 0.4036286627736071;
    s0 = 0.3901759381306066;
    s1 = 0.5501596788096901;
    s2 = 0.4524276954209931;
    F_wkspc<MP0, 2> w;
    set_args<MP0, 3, P110, P101, P011>(w, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor_bv<MP0, 3, P110, P101, P011> func {w};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], -2.775557561562891e-17, 2.22045e-16);
    ASSERT_NEAR(info.lambda[1], -8.881784197001252e-16, 9e-16);
  }
}

TEST (numopt, sqp_bary_works_with_mp1_111) {
  double u0, u1, u2, h, s, s0, s1, s2;
  bool error;
  {
    u0 = 0.4248584117046258;
    u1 = 0.4603291735753289;
    u2 = 0.5721682403965839;
    h = 0.2975553841511184;
    s = 0.7064072275375605;
    s0 = 0.2435733726809508;
    s1 = 0.477175202496374;
    s2 = 0.2656191251853977;
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3, P100, P010, P001>(w, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor_bv<MP1, 3, P100, P010, P001> func {w};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0.3535224569679952, 2.22045e-16);
    ASSERT_NEAR(info.lambda[1], 0, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_mp1_123) {
  double u0, u1, u2, h, s, s0, s1, s2;
  bool error;
  {
    u0 = 0.9463249898054826;
    u1 = 0.9762491369597006;
    u2 = 0.9682220862053195;
    h = 0.03918448664758301;
    s = 0.1838429444657746;
    s0 = 0.4979488150189469;
    s1 = 0.5182403290268277;
    s2 = 0.5369077169939814;
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3, P100, P110, P111>(w, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor_bv<MP1, 3, P100, P110, P111> func {w};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], -8.881784197001252e-16, 9e-16);
    ASSERT_NEAR(info.lambda[1], -7.771561172376096e-16, 8e-16);
  }
}

TEST (numopt, sqp_bary_works_with_mp1_222) {
  double u0, u1, u2, h, s, s0, s1, s2;
  bool error;
  {
    u0 = 0.4634892477622432;
    u1 = 0.4643570680407132;
    u2 = 0.5485765253534179;
    h = 0.09298892687067795;
    s = 0.642741739133104;
    s0 = 0.001419058202402956;
    s1 = 0.004244552084378484;
    s2 = 0.02080448055205524;
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3, P110, P101, P011>(w, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor_bv<MP1, 3, P110, P101, P011> func {w};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0.4790032956983284, 7e-16);
    ASSERT_NEAR(info.lambda[1], 0, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_rhr_111) {
  double u0, u1, u2, h, s, s0, s1, s2;
  bool error;
  {
    u0 = 0.4365547823722682;
    u1 = 0.4580440444309234;
    u2 = 0.4582268818076987;
    h = 0.4366566416352912;
    s = 0.09110017541305393;
    s0 = 0.5940370314441212;
    s1 = 0.6993079853260437;
    s2 = 0.9614264378869783;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3, P100, P010, P001>(w, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor_bv<RHR, 3, P100, P010, P001> func {w};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0.2183590041368041, 2e-13);
    ASSERT_NEAR(info.lambda[1], 0.2153991075243049, 2e-13);
  }
}

TEST (numopt, sqp_bary_works_with_rhr_123) {
  double u0, u1, u2, h, s, s0, s1, s2;
  bool error;
  {
    u0 = 0.9046664796901251;
    u1 = 0.9844082946380126;
    u2 = 0.9442484229394705;
    h = 0.09016600216896353;
    s = 0.7817226129171663;
    s0 = 0.1484650229434085;
    s1 = 0.2043513471851703;
    s2 = 0.171964418184969;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3, P100, P110, P111>(w, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor_bv<RHR, 3, P100, P110, P111> func {w};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0, 2.22045e-16);
    ASSERT_NEAR(info.lambda[1], 0, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_rhr_222) {
  double u0, u1, u2, h, s, s0, s1, s2;
  bool error;
  {
    u0 = 0.6319307979774571;
    u1 = 1.027881046467556;
    u2 = 0.8567754790043727;
    h = 0.4018833980080856;
    s = 0.933591915842209;
    s0 = 0.7203432060007481;
    s1 = 0.9148702499040858;
    s2 = 0.9771591929603007;
    F_wkspc<RHR, 2> w;
    set_args<RHR, 3, P110, P101, P011>(w, u0, u1, u2, s, s0, s1, s2, h);
    cost_functor_bv<RHR, 3, P110, P101, P011> func {w};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 5.551115123125783e-17, 2.22045e-16);
    ASSERT_NEAR(info.lambda[1], 0.09484711766285059, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_mp0_fac) {
  double u0, u1, u2, h, s, s0, s1, s2, sf, p0[3], p1[3], p2[3], pf[3];
  bool error;
  {
    u0 = 0.5201903182114772;
    u1 = 1.423276387824696;
    u2 = 0.5898718928250356;
    h = 0.9468166671615945;
    s = 0.207031948073074;
    s0 = 0.7750278139819674;
    s1 = 1.640596080051717;
    s2 = 1.51595981013218;
    sf = 0.295534196164876;
    p0[0] = 0.1518457223144846;
    p0[1] = 0.8479105223137822;
    p0[2] = 0.7848545910931496;
    p1[0] = 0.2708315021513173;
    p1[1] = 0.2278107048161337;
    p1[2] = 0.3210232170889511;
    p2[0] = 0.8295618043991675;
    p2[1] = 0.8221821946009759;
    p2[2] = 0.5706828504818783;
    pf[0] = 0.2404774846793727;
    pf[1] = -0.8215003988470204;
    pf[2] = 0.9931121262131963;
    F_fac_wkspc<MP0, 2> w;
    set_args<MP0, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, pf, sf);
    cost_functor_fac<MP0, 3> func {w, p0, p1, p2, pf};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0, 2.22045e-16);
    ASSERT_NEAR(info.lambda[1], 0.06325936565896506, 4e-16);
  }
  {
    u0 = 0;
    u1 = 1;
    u2 = sqrt2;
    s = 1;
    s0 = 1;
    s1 = 1;
    s2 = 1;
    sf = 1;
    h = 1;
    p0[0] = 1;
    p0[1] = 0;
    p0[2] = -1;
    p1[0] = 1;
    p1[1] = 0;
    p1[2] = 0;
    p2[0] = 1;
    p2[1] = 1;
    p2[2] = 0;
    pf[0] = 1;
    pf[1] = 0;
    pf[2] = -1;
    F_fac_wkspc<MP0, 2> w;
    set_args<MP0, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, pf, sf);
    cost_functor_fac<MP0, 3> func {w, p0, p1, p2, pf};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0, 2.22045e-16);
    ASSERT_NEAR(info.lambda[1], 0, 2.22045e-16);
    ASSERT_NEAR(info.value, sqrt2, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_mp1_fac) {
  double u0, u1, u2, h, s, s0, s1, s2, sf, p0[3], p1[3], p2[3], pf[3];
  bool error;
  {
    u0 = 0.8901234779616869;
    u1 = 1.242794667059232;
    u2 = 1.074934814085041;
    h = 0.6548828984496685;
    s = 0.9759575178799641;
    s0 = 0.03642551552548012;
    s1 = 0.2500775071346921;
    s2 = 0.673635497771426;
    sf = 0.3650326253057474;
    p0[0] = 0.3091496188959173;
    p0[1] = 0.1209123845806283;
    p0[2] = 0.9157657040407896;
    p1[0] = 0.1354782056815019;
    p1[1] = 0.3321178918283407;
    p1[2] = 0.8974798923998842;
    p2[0] = 0.4996487778472734;
    p2[1] = 0.6152882422150567;
    p2[2] = 0.5831329680886482;
    pf[0] = 0.6808677941875395;
    pf[1] = -1.287415913860736;
    pf[2] = 0.2180848525247276;
    F_fac_wkspc<MP1, 2> w;
    set_args<MP1, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, pf, sf);
    cost_functor_fac<MP1, 3> func {w, p0, p1, p2, pf};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0, 2.22045e-16);
    ASSERT_NEAR(info.lambda[1], 0, 2.22045e-16);
  }
}

TEST (numopt, sqp_bary_works_with_rhr_fac) {
  double u0, u1, u2, h, s, s0, s1, s2, sf, p0[3], p1[3], p2[3], pf[3];
  bool error;
  {
    u0 = 0.2297011978711195;
    u1 = 0.2673272949893583;
    u2 = 0.3323685751020444;
    h = 0.3302022425140205;
    s = 0.2284323221357374;
    s0 = 0.6519971672170032;
    s1 = 0.6738433887243777;
    s2 = 0.742945226505809;
    sf = 0.2818202375608668;
    p0[0] = 0.8800662597260562;
    p0[1] = 0.4443303586992347;
    p0[2] = 0.7559141209079612;
    p1[0] = 0.6032963758153485;
    p1[1] = 0.7832659374040363;
    p1[2] = 0.1139306386784259;
    p2[0] = 0.9785638851599779;
    p2[1] = 0.848596675479549;
    p2[2] = 0.05064649767054341;
    pf[0] = -0.1764880040665571;
    pf[1] = -0.6520821914757466;
    pf[2] = 0.1565338036926401;
    F_fac_wkspc<RHR, 2> w;
    set_args<RHR, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, pf, sf);
    cost_functor_fac<RHR, 3> func {w, p0, p1, p2, pf};
    updates::info<2> info;
    sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);
    ASSERT_FALSE(error);
    ASSERT_NEAR(info.lambda[0], 0.2424829348139743, 2e-15);
    ASSERT_NEAR(info.lambda[1], 0, 2.22045e-16);
  }
}
