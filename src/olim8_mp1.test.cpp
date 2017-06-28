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

BOOST_AUTO_TEST_CASE (origin_test) {
  size_t M = 5, N = 5;
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

BOOST_AUTO_TEST_CASE (sf1_single_row_test) {
  size_t N = 1001;
  double h = 1.0/N;
  olim8_mp1 m {1, N, h, sf1};
  m.add_boundary_node(0, 0);
  m.run();
  for (size_t j = N - 10; j < N; ++j) {
    double U = m.get_value(0, j);
    double u = sf1_soln(h*j, 0);
    printf("U = %g, u = %g (err: %g)\n", U, u, fabs(u - U)/fabs(u));
    // BOOST_TEST(u == U, boost::test_tools::tolerance(1e-2));
  }
}

// BOOST_AUTO_TEST_CASE (sf1_test) {
//   size_t M = 101, N = M;
//   double h = 1.0/(M - 1);
//   olim8_mp1 m {M, N, h, sf1};
//   m.add_boundary_node(0, 0);
//   m.run();
//   for (size_t i = M - 3; i < M; ++i) {
//     for (size_t j = N - 3; j < N; ++j) {
//       double U = m.get_value(i, j);
//       double u = sf1_soln(h*j, h*i);
//       printf("%g\n", fabs(u - U)/fabs(u));
//       BOOST_TEST(u == U, boost::test_tools::tolerance(1e-1));
//     }
//   }
// }

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
