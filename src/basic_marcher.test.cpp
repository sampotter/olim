#define BOOST_TEST_MODULE basic_marcher

#include <boost/test/included/unit_test.hpp>

#include "basic_marcher.hpp"
#include "speed_funcs.hpp"

BOOST_AUTO_TEST_CASE (trivial_case_works) {
  basic_marcher m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
}

BOOST_AUTO_TEST_CASE (adjacent_update_works) {
  basic_marcher m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 0.5);
}

BOOST_AUTO_TEST_CASE (neighboring_values_are_correct) {
  basic_marcher m {3, 3, 1};
  m.add_boundary_node(1, 1);
  m.run();
  double gt[] = {
    1.707106781186547,
    1,
    1.707106781186547,
    1,
    0,
    1,
    1.707106781186547,
    1,
    1.707106781186547
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

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
