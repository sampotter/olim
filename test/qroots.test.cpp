#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "qroots.hpp"
#include "speed_funcs.hpp"

TEST (qroots, sigma_works_1) {
  double a[] = {
    -0.0019249999999999998,
    -0.0010000000000000002,
    0.00087500000000000078,
    -0.0020000000000000005,
    0.00040000000000000018
  };
  double b[] = {
    -0.0010000000000000002,
    0.0017500000000000016,
    -0.0060000000000000019,
    0.0016000000000000007
  };
  double c[] = {
    0.0022374999999999999,
    0.00020312499999999986,
    0.0014374999999999998
  };
  double d[] = {
    -0.0086910396975425352,
    -0.0051202079395085099
  };
  double e[] = {
    -0.006034391692523379
  };
  double * polys[] = {a, b, c, d, e};

  ASSERT_TRUE(sigma(polys, 0.1) == 2);
}

TEST (qroots, sigma_works_2) {
  double a[] = {
    -0.99000016666583334,
    0.0098791278581576019,
    -0.013920963747804107,
    0.019758255716315204,
    0.000099996666711110551
  };
  double b[] = {
    0.0098791278581576019,
    -0.027841927495608213,
    0.059274767148945612,
    0.0003999866668444422
  };
  double c[] = {
    1.1120006923430827,
    -0.35123825863789943,
    0.738963635937399
  };
  double d[] = {
    0.079604292893865752,
    -0.00102434857077549
  };
  double e[] = {
    -4436.5455201568111
  };
  double * polys[] = {a, b, c, d, e};
  int lhs = sigma(polys, 0), oldlhs = oldsigma(polys, 0);
  int rhs = sigma(polys, 1), oldrhs = oldsigma(polys, 1);
  ASSERT_TRUE(lhs == 2);
  ASSERT_TRUE(oldlhs == 2);
  ASSERT_TRUE(rhs == 2);
  ASSERT_TRUE(oldrhs == 2);
}

TEST (qroots, sigma_works_3) {
  double a[] = {
    0,
    0.000050000000000000002,
    0.0001,
    0.0001,
    0.0001
  };
  double b[] = {
    0.000050000000000000002,
    0.00020000000000000001,
    0.00030000000000000003,
    0.00040000000000000002
  };
  double c[] = {
    0.0000031249999999999997,
    -0.000025000000000000005,
    -0.000031250000000000001
  };
  double d[] = {
    -0.000048000000000000001,
    -0.00017600000000000005
  };
  double e[] = {
    -0.0000076188016528925621
  };
  double * polys[] = {a, b, c, d, e};
  double eps = std::numeric_limits<double>::epsilon();
  int lhs = sigma(polys, 0 + eps), oldlhs = oldsigma(polys, 1);
  int rhs = sigma(polys, 0 + eps), oldrhs = oldsigma(polys, 1);
  ASSERT_TRUE(lhs == oldlhs);
  ASSERT_TRUE(rhs == oldrhs);
}

TEST (qroots, sturm_works_1) {
  double a[] = {
    -0.0019249999999999998,
    -0.0010000000000000002,
    0.00087500000000000078,
    -0.0020000000000000005,
    0.00040000000000000018
  };
  double b[] = {
    -0.0010000000000000002,
    0.0017500000000000016,
    -0.0060000000000000019,
    0.0016000000000000007
  };
  double c[] = {
    0.0022374999999999999,
    0.00020312499999999986,
    0.0014374999999999998
  };
  double d[] = {
    -0.0086910396975425352,
    -0.0051202079395085099
  };
  double e[] = {
    -0.006034391692523379
  };
  double * polys[] = {a, b, c, d, e};
  ASSERT_TRUE(sturm(polys, 0, 1) == 0);
}

TEST (qroots, sturm_works_2) {
  double a[] = {
    0.0044039823074159008,
    -0.018607964614831808,
    0.018607964614831819,
    -0.000000000000000015543122344752199,
    0.0000000000000000000000000000000030814879110195774
  };
  double b[] = {
    -0.018607964614831808,
    0.037215929229663637,
    -0.000000000000000046629367034256598,
    0.000000000000000000000000000000012325951644078309
  };
  double c[] = {
    5866197575383.5918,
    -11732395150767.186,
    0.0053960176925841059
  };
  double d[] = {
    -0.0029492948003932741,
    0.00589858960078653
  };
  double e[] = {
    0.0185546875
  };
  double * polys[] = {a, b, c, d, e};
  ASSERT_TRUE(sturm(polys, 0, 1) == 2);
  ASSERT_TRUE(sturm(polys, 0, 0.5) == 1);
  ASSERT_TRUE(sturm(polys, 0.5, 1) == 1);
  ASSERT_TRUE(sturm(polys, 0.4, 0.5) == 0);
  ASSERT_TRUE(sturm(polys, 0.49999999999999994, 0.5) == 0);
}

TEST (qroots, sturm_works_3) {
  double a[] = {
    0.000075000000000000007,
    -0.00020000000000000001,
    0.000475,
    -0.00040000000000000002,
    0.00040000000000000002
  };
  double b[] = {
    -0.00020000000000000001,
    0.00094999999999999999,
    -0.0012000000000000001,
    0.0016000000000000001
  };
  double c[] = {
    -0.000062500000000000001,
    0.000090625000000000021,
    -0.00016249999999999999
  };
  double d[] = {
    0.000081656804733727889,
    -0.0013937869822485207
  };
  double e[] = {
    0.000057748375077573313
  };
  double * polys[] = {a, b, c, d, e};
  ASSERT_TRUE(sturm(polys, 0, 1) == 0);
}

TEST (qroots, sturm_works_4) {
  double a[] = {
    -0.99476271779816383,
    -0.0040797065327272484,
    0.80278768800205891,
    -0.0081594130654544968,
    0.00001715476311934851
  };
  double b[] = {
    -0.0040797065327272484,
    1.6055753760041178,
    -0.02447823919636349,
    0.000068619052477394041
  };
  double c[] = {
    0.2887299652279402,
    -47.72622322373811,
    0.32627601410407686
  };
  double d[] = {
    -0.0086994542995620129,
    0.50672200724493743
  };
  double e[] = {
    0.5305424594038346
  };
  double * polys[] = {a, b, c, d, e};
  ASSERT_TRUE(sturm(polys, 0, 1) == 0);
}

// Some helper functions for testing the roots returned by qroots:

static int num_roots(double const * roots) {
  int i = 0;
  while (roots[i] != -1) i++;
  return i;
}

static bool contains_root(double const * roots, double root,
                          double tol = 1e-10) {
  return fabs(roots[0] - root)/fabs(root) < tol ||
    fabs(roots[1] - root)/fabs(root) < tol ||
    fabs(roots[2] - root)/fabs(root) < tol ||
    fabs(roots[3] - root)/fabs(root) < tol;
}

TEST (qroots, qroots_works_1) {
  double a[5] = {-2.25, 11.8125, -10.8125, -2, 1};
  double roots[4] = {-1, -1, -1, -1};
  qroots(a, roots, 0, 1);
  ASSERT_TRUE(num_roots(roots) == 2);
  ASSERT_TRUE(contains_root(roots, 0.25));
  ASSERT_TRUE(contains_root(roots, 0.75));
}

TEST (qroots, qroots_works_2) {
  double a[5] = { -0.803056, 2.86611, -1.10611, -2.56, -0.64};
  double roots[4] = {-1, -1, -1, -1};
  qroots(a, roots, 0, 1);
  ASSERT_TRUE(num_roots(roots) == 2);
  ASSERT_TRUE(contains_root(roots, 0.432017, 1e-6));
  ASSERT_TRUE(contains_root(roots, 0.4828, 1e-5));
}
