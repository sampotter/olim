#define BOOST_TEST_MODULE basic_marcher

#include <boost/test/included/unit_test.hpp>

#include "basic_marcher.hpp"

BOOST_AUTO_TEST_CASE (trivial_case_works) {
  basic_marcher m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
};

BOOST_AUTO_TEST_CASE (adjacent_update_works) {
  basic_marcher m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 1), 0.5);
};

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
