#define BOOST_TEST_MODULE olim8_rhr

#include <boost/test/included/unit_test.hpp>

#include "olim8_rhr.hpp"

#include "speed_funcs.hpp"

BOOST_AUTO_TEST_CASE (trivial_case_works) {
  olim8_rhr m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
};

BOOST_AUTO_TEST_CASE (adjacent_update_works) {
  olim8_rhr m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 0.5);
};

BOOST_AUTO_TEST_CASE (neighboring_values_are_correct) {
  olim8_rhr m {3, 3, 1};
  m.add_boundary_node(1, 1);
  m.run();
  double gt[] = {
    std::sqrt(2),
    1,
    std::sqrt(2),
    1,
    0,
    1,
    std::sqrt(2),
    1,
    std::sqrt(2)
  };
  int k = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      BOOST_TEST(
        gt[k++] == m.get_value(i, j),
        boost::test_tools::tolerance(1e-15));
    }
  }
}

BOOST_AUTO_TEST_CASE (slightly_more_involved) {
  olim8_rhr m {2, 2, 1};
  m.add_boundary_node(0, 0);
  m.add_boundary_node(1, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(0, 1), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 1), 1.0);
}

BOOST_AUTO_TEST_CASE (slightly_more_involved_2) {
  olim8_rhr m {2, 2, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(0, 1), 1.0);
}

BOOST_AUTO_TEST_CASE (run_on_full_neighborhood) {
  olim8_rhr m {3, 3, 1};
  m.add_boundary_node(0, 0);
  m.add_boundary_node(0, 1);
  m.add_boundary_node(0, 2);
  m.add_boundary_node(1, 0);
  m.add_boundary_node(1, 2);
  m.add_boundary_node(2, 0);
  m.add_boundary_node(2, 1);
  m.add_boundary_node(2, 2);
  m.run();
}

BOOST_AUTO_TEST_CASE (maria_test, *boost::unit_test::tolerance(1e-15)) {
  olim8_rhr m {3, 2, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(0, 1), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 1), std::sqrt(2.0));
  BOOST_CHECK_EQUAL(m.get_value(2, 0), 2.0);
  BOOST_CHECK_EQUAL(m.get_value(2, 1), 2.3243932834975496);
}

BOOST_AUTO_TEST_CASE (origin_test) {
  int M = 5, N = 5;
  double h = 0.1, x0 = h*(N - 1)/2., y0 = h*(M - 1)/2.;
  olim8_rhr m {M, N, h, default_speed_func, x0, y0};
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
  olim8_rhr m {1, N, h, s1};
  m.add_boundary_node(0, 0);
  m.run();
  double maxerr = 0;
  for (int j = 0; j < N; ++j) {
    double U = m.get_value(0, j);
    double u = f1(h*j, 0);
    maxerr = std::max(maxerr, std::fabs(U - u)/std::fabs(u));
  }
  std::cout << maxerr << std::endl;
}

BOOST_AUTO_TEST_CASE (s1_test1) {
  int M = 1001, N = M;
  double h = 1.0/(M - 1);
  olim8_rhr m {M, N, h, s1};
  m.add_boundary_node(0, 0);
  m.run();
  double maxerr = 0;
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f1(h*j, h*i);
      if (u != 0) {
        maxerr = std::max(maxerr, std::fabs(U - u)/std::fabs(u));
      }
    }
  }
  std::cout << maxerr << std::endl;
}

BOOST_AUTO_TEST_CASE (s1_test2) {
  int M = 201, N = M;
  double h = 2.0/(M - 1);
  olim8_rhr m {M, N, h, s1};
  m.add_boundary_node(100, 100);
  m.run();
  double maxerr = 0;
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f1(h*j, h*i);
      if (u != 0) {
        maxerr = std::max(maxerr, std::fabs(U - u)/std::fabs(u));
      }
    }
  }
  std::cout << maxerr << std::endl;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
