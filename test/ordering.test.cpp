#include <gtest/gtest.h>

#include "olim.hpp"

TEST (ordering, row_major_works) {
  using olim_t = olim4_mp0<ordering::ROW_MAJOR>;
  typename olim_t::ivec dims {2, 3};
  double h = 1;
  double s[6] = {1, 2, 3, 4, 5, 6};
  olim_t olim {dims, h, s};
  olim.add_src({1, 0});
  olim.add_bd({0, 2});
  ASSERT_TRUE(isinf(olim.get_U({0, 0})));
  ASSERT_TRUE(isinf(olim.get_U({0, 1})));
  ASSERT_TRUE(isinf(olim.get_U({0, 2})));
  ASSERT_EQ(olim.get_U({1, 0}), 0);
  ASSERT_TRUE(isinf(olim.get_U({1, 1})));
  ASSERT_TRUE(isinf(olim.get_U({1, 2})));
}
