#define BOOST_TEST_MODULE olim_8pt_rhr

#include <boost/test/included/unit_test.hpp>

#include "olim_8pt_rhr.hpp"

BOOST_AUTO_TEST_CASE (trivial_case_works) {
  olim_8pt_rhr m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
};

BOOST_AUTO_TEST_CASE (adjacent_update_works) {
  olim_8pt_rhr m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 0.5);
};

BOOST_AUTO_TEST_CASE (neighboring_values_are_correct) {
  olim_8pt_rhr m {3, 3, 1};
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
  olim_8pt_rhr m {2, 2, 1};
  m.add_boundary_node(0, 0);
  m.add_boundary_node(1, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(0, 1), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 1), 1.0);
}

BOOST_AUTO_TEST_CASE (slightly_more_involved_2) {
  olim_8pt_rhr m {2, 2, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(0, 1), 1.0);
}

BOOST_AUTO_TEST_CASE (run_on_full_neighborhood) {
  olim_8pt_rhr m {3, 3, 1};
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
  olim_8pt_rhr m {3, 2, 1};
  m.add_boundary_node(0, 0);
  m.run();
  BOOST_CHECK_EQUAL(m.get_value(0, 0), 0.0);
  BOOST_CHECK_EQUAL(m.get_value(0, 1), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 0), 1.0);
  BOOST_CHECK_EQUAL(m.get_value(1, 1), std::sqrt(2.0));
  BOOST_CHECK_EQUAL(m.get_value(2, 0), 2.0);
  BOOST_CHECK_EQUAL(m.get_value(2, 1), 2.3243932834975496);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
