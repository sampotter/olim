#define BOOST_TEST_MODULE solim4_mp0

#include <boost/test/included/unit_test.hpp>

#include <cstdio>

#include "solim4_mp0.hpp"

BOOST_AUTO_TEST_CASE (trying_it_out) {
  solim4_mp0 m {11, 11};
  m.add_boundary_node(5, 5);
  m.run();
  for (int i = 0; i <= 10; ++i) {
    for (int j = 0; j <= 10; ++j) {
      double value = m.get_value(i, j), lambda = m.get_lambda(i, j);
      printf("%d, %d: %g, %g\n", i, j, value, lambda);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
