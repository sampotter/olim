#define BOOST_TEST_MODULE olim8_mp1

#include <boost/test/included/unit_test.hpp>

#include "olim8_mp1.hpp"
#include "speed_funcs.hpp"

BOOST_AUTO_TEST_CASE (trivial_case_works) {
  olim8_mp1 m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
}

BOOST_AUTO_TEST_CASE (adjacent_update_works) {
  olim8_mp1 m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 0.5);
}

BOOST_AUTO_TEST_CASE (neighboring_values_are_correct) {
  olim8_mp1 m {3, 3, 1};
  m.add_boundary_node(1, 1);
  m.run();
  double root2 = std::sqrt(2);
  double gt[] = {root2, 1, root2, 1, 0, 1, root2, 1, root2};
  int k = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      BOOST_TEST(
        gt[k++] == m.get_value(i, j),
        boost::test_tools::tolerance(1e-15));
    }
  }
}

BOOST_AUTO_TEST_CASE (origin_test) {
  int M = 5, N = 5;
  double h = 0.1, x0 = h*(N - 1)/2., y0 = h*(M - 1)/2.;
  olim8_mp1 m {M, N, h, default_speed_func, x0, y0};
  m.add_boundary_node(2, 2);
  m.run();
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      BOOST_TEST(
        m.get_value(i, j) == default_speed_func_soln(h*j - x0, h*i - y0),
        boost::test_tools::tolerance(4e-2));
    }
  }
}

BOOST_AUTO_TEST_CASE (s1_single_row_test) {
  int N = 1001;
  double h = 1.0/(N - 1);
  olim8_mp1 m {1, N, h, s1};
  m.add_boundary_node(0, 0);
  m.run();
  for (int j = N - 10; j < N; ++j) {
    double U = m.get_value(0, j);
    double u = f1(h*j, 0);
    BOOST_TEST(u == U, boost::test_tools::tolerance(1e-2));
  }
}

BOOST_AUTO_TEST_CASE (s1_test) {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp1 m {M, N, h, s1};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f1(h*j, h*i);
      BOOST_TEST(u == U, boost::test_tools::tolerance(1e-1));
    }
  }
}

BOOST_AUTO_TEST_CASE (s2_test) {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp1 m {M, N, h, s2};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f2(h*j, h*i);
      BOOST_TEST(u == U, boost::test_tools::tolerance(1e-1));
    }
  }
}

BOOST_AUTO_TEST_CASE (s3_test) {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp1 m {M, N, h, s3};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f3(h*j, h*i);
      BOOST_TEST(u == U, boost::test_tools::tolerance(1e-1));
    }
  }
}

BOOST_AUTO_TEST_CASE (s4_test) {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp1 m {M, N, h, s4};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f4(h*j, h*i);
      BOOST_TEST(u == U, boost::test_tools::tolerance(1e-1));
    }
  }
}

BOOST_AUTO_TEST_CASE (s5_test) {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp1 m {M, N, h, s5};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f5(h*j, h*i);
      BOOST_TEST(u == U, boost::test_tools::tolerance(1e-1));
    }
  }
}

BOOST_AUTO_TEST_CASE (s6_test) {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp1 m {M, N, h, s6};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f6(h*j, h*i);
      BOOST_TEST(u == U, boost::test_tools::tolerance(1e-1));
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
