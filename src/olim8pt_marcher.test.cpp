#define BOOST_TEST_MODULE olim8pt_marcher

#include <boost/test/included/unit_test.hpp>

#include "olim8pt_marcher.hpp"

BOOST_AUTO_TEST_CASE (trivial_case_works) {
  olim8pt_marcher m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
};

BOOST_AUTO_TEST_CASE (adjacent_update_works) {
  olim8pt_marcher m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 0.5);
};

BOOST_AUTO_TEST_CASE (slightly_more_involved) {
  olim8pt_marcher m {2, 2, 1};
  m.add_boundary_node(0, 0);
  m.add_boundary_node(1, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(0, 1), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 1), 1.0);
}

BOOST_AUTO_TEST_CASE (slightly_more_involved_2) {
  olim8pt_marcher m {2, 2, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(0, 1), 1.0);
}

BOOST_AUTO_TEST_CASE (run_on_full_neighborhood) {
  olim8pt_marcher m {3, 3, 1};
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

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
