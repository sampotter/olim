#include <gtest/gtest.h>

#include "olim.hpp"

TEST (marcher, fractional_add_boundaries_is_correct_with_constant_slowness) {
  double S[4] = {1, 1, 1, 1};
  olim4_rhr o {2, 2, 1, S};
  o.add_boundary_node(3./8, 1./4, 1.);
  ASSERT_DOUBLE_EQ(o.get_value(0, 0), 0.45069390943299864);
  ASSERT_DOUBLE_EQ(o.get_value(0, 1), 0.673145600891813);
  ASSERT_DOUBLE_EQ(o.get_value(1, 0), 0.8385254915624212);
  ASSERT_DOUBLE_EQ(o.get_value(1, 1), 0.9762812094883317);
}

TEST (marcher, fractional_add_boundaries_is_correct_with_nonconstant_slowness) {
  double s = 0.9, S[4] = {1.1, 1.2, 1.3, 1.4}, x = 3./16, y = 1./8, h = 0.5;
  double l[4] = {
    0.45069390943299864, 0.673145600891813,
    0.8385254915624212, 0.9762812094883317
  };
  {
    olim4_rhr o {2, 2, h, S};
    o.add_boundary_node(x, y, s);
    ASSERT_DOUBLE_EQ(o.get_value(0, 0), S[0]*h*l[0]);
    ASSERT_DOUBLE_EQ(o.get_value(0, 1), S[1]*h*l[1]);
    ASSERT_DOUBLE_EQ(o.get_value(1, 0), S[2]*h*l[2]);
    ASSERT_DOUBLE_EQ(o.get_value(1, 1), S[3]*h*l[3]);
  }
  {
    olim4_mp0 o {2, 2, h, S};
    o.add_boundary_node(x, y, s);
    ASSERT_DOUBLE_EQ(o.get_value(0, 0), (s + S[0])*h*l[0]/2);
    ASSERT_DOUBLE_EQ(o.get_value(0, 1), (s + S[1])*h*l[1]/2);
    ASSERT_DOUBLE_EQ(o.get_value(1, 0), (s + S[2])*h*l[2]/2);
    ASSERT_DOUBLE_EQ(o.get_value(1, 1), (s + S[3])*h*l[3]/2);
  }
}
