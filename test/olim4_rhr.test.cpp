#include "basic_marcher.hpp"
#include "olim.hpp"
#include "olim.test.common.hpp"

using olim_t = olim4_rhr;

TEST (olim4_rhr, trivial_case_works) {
  ASSERT_TRUE(trivial_case_works<olim_t>());
}

TEST (olim4_rhr, adjacent_update_works) {
  ASSERT_TRUE(adjacent_update_works<olim_t>());
}

TEST (olim4_rhr, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim_t>(1.0 + 1.0/sqrt(2)));
}

TEST (olim4_rhr, correct_corners_in_limit) {
  ASSERT_TRUE(correct_corners_in_limit<olim_t>(101, 2.7e-2));
}

TEST (olim4_rhr, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) default_speed_func));
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) s1));
}

TEST (olim4_rhr, agrees_with_basic_marcher) {
  for (auto & s: speed_funcs) {
    auto res = olims_agree<olim_t, basic_marcher>(s);
    ASSERT_TRUE(res);
  }
}

TEST (olim4_rhr, solution_is_exact_in_factored_region) {
  olim_t o {3, 3, 1, (speed_func) default_speed_func, 1., 1.};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == 1 && j == 1) continue;
      o.set_node_parent(i, j, 1, 1);
    }
  }
  o.add_boundary_node(1, 1);
  o.run();

  for (int i = 0; i < 3; ++i) {
    double y = i - 1;
    for (int j = 0; j < 3; ++j) {
      double x = j - 1;
      double U = o.get_value(i, j);
      ASSERT_EQ(U, std::hypot(x, y));
    }
  }
}
