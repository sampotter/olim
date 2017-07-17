#define BOOST_TEST_MODULE solim4_mp0

#include <boost/test/included/unit_test.hpp>

#include <cstdio>

#include "solim4_mp0.hpp"

BOOST_AUTO_TEST_CASE (trying_it_out) {
  int M = 31, N = 31;
  solim4_mp0 m {M, N};
  m.add_boundary_node(5, 5);
  m.add_boundary_node(27, 3);
  m.add_boundary_node(11, 19);
  m.run();
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double value = m.get_value(i, j), lambda = m.get_lambda(i, j);
      auto const parents = m.get_parent_nodes(i, j);
      auto const n0 = static_cast<node *>(parents.first);
      auto const n1 = static_cast<node *>(parents.second);
      if (n0 && n1) {
        printf(
          "%d, %d: T = %g, lam = %g, i0 = %d, j0 = %d, i1 = %d, j1 = %d\n",
          i, j, value, lambda, n0->get_i(), n0->get_j(), n1->get_i(),
          n1->get_j());
      } else if (n0) {
        printf("%d, %d: T = %g, i0 = %d, j0 = %d\n", i, j, value, n0->get_i(),
               n0->get_j());
      } else {
        printf("%d, %d: T = %g\n", i, j, value);
      }
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
