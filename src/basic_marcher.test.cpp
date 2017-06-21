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

BOOST_AUTO_TEST_CASE (sf1_values_are_correct) {
  basic_marcher m {5, 5, 0.1, sf1};
  m.add_boundary_node(2, 2);
  m.run();
  double gt[] = {
    0.4551098033382107,
    0.3763043590877354,
    0.3115278386858286,
    0.4017273107926478,
    0.5164623464645542,
    0.3397867999746174,
    0.2396828846632308,
    0.1495134660269984,
    0.2571280920182202,
    0.3850073790874829,
    0.2551369471830275,
    0.1316278775317397,
    0,
    0.144029241953726,
    0.2902839752847525,
    0.3169807447134633,
    0.217968740635256,
    0.1295135310607615,
    0.2309905205398307,
    0.3515496483603887,
    0.3942743416649593,
    0.3151037675674407,
    0.2508406384704259,
    0.3284435944832341,
    0.4280955167506931
  };
  int k = 0;
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
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
