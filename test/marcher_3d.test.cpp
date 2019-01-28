#include <gtest/gtest.h>

#include "olim3d.hpp"

double l3(double x, double y, double z) {
  return sqrt(x*x + y*y + z*z);
}

TEST (marcher_3d, fractional_add_boundaries_is_correct_with_constant_slowness) {
  double S[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  double x = 2./8, y = 1./8, z = 3./8;
  double l[8] = {
    l3(x,     y,     z),     // 0, 0, 0
    l3(x,     y,     1 - z), // 0, 0, 1
    l3(x,     1 - y, z),     // 0, 1, 0
    l3(1 - x, y,     z),     // 1, 0, 0
    l3(x,     1 - y, 1 - z), // 0, 1, 1
    l3(1 - x, y,     1 - z), // 1, 0, 1
    l3(1 - x, 1 - y, z),     // 1, 1, 0
    l3(1 - x, 1 - y, 1 - z)  // 1, 1, 1
  };
  olim6_rhr o {{2, 2, 2}, 1, S};
  o.add_boundary_node(vec3<double> {x, y, z}, 1.);
  ASSERT_DOUBLE_EQ(o.get_value({0, 0, 0}), l[0]);
  ASSERT_DOUBLE_EQ(o.get_value({0, 0, 1}), l[1]);
  ASSERT_DOUBLE_EQ(o.get_value({0, 1, 0}), l[2]);
  ASSERT_DOUBLE_EQ(o.get_value({1, 0, 0}), l[3]);
  ASSERT_DOUBLE_EQ(o.get_value({0, 1, 1}), l[4]);
  ASSERT_DOUBLE_EQ(o.get_value({1, 0, 1}), l[5]);
  ASSERT_DOUBLE_EQ(o.get_value({1, 1, 0}), l[6]);
  ASSERT_DOUBLE_EQ(o.get_value({1, 1, 1}), l[7]);
}

TEST (marcher_3d, fractional_add_boundaries_is_correct_with_nonconstant_slowness) {
  double s = 0.9, S[8] = {1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35}, h = 0.5;
  double x = 2./8, y = 1./8, z = 3./8;
  double l[8] = {
    l3(x/h,     y/h,     z/h),     // 0, 0, 0
    l3(1 - x/h, y/h,     z/h),     // 1, 0, 0
    l3(x/h,     1 - y/h, z/h),     // 0, 1, 0
    l3(1 - x/h, 1 - y/h, z/h),     // 1, 1, 0
    l3(x/h,     y/h,     1 - z/h), // 0, 0, 1
    l3(1 - x/h, y/h,     1 - z/h), // 1, 0, 1
    l3(x/h,     1 - y/h, 1 - z/h), // 0, 1, 1
    l3(1 - x/h, 1 - y/h, 1 - z/h)  // 1, 1, 1
  };
  {
    olim6_rhr o {{2, 2, 2}, h, S};
    o.add_boundary_node(vec3<double> {x, y, z}, s);
    ASSERT_DOUBLE_EQ(o.get_value({0, 0, 0}), S[0]*h*l[0]);
    ASSERT_DOUBLE_EQ(o.get_value({1, 0, 0}), S[1]*h*l[1]);
    ASSERT_DOUBLE_EQ(o.get_value({0, 1, 0}), S[2]*h*l[2]);
    ASSERT_DOUBLE_EQ(o.get_value({1, 1, 0}), S[3]*h*l[3]);
    ASSERT_DOUBLE_EQ(o.get_value({0, 0, 1}), S[4]*h*l[4]);
    ASSERT_DOUBLE_EQ(o.get_value({1, 0, 1}), S[5]*h*l[5]);
    ASSERT_DOUBLE_EQ(o.get_value({0, 1, 1}), S[6]*h*l[6]);
    ASSERT_DOUBLE_EQ(o.get_value({1, 1, 1}), S[7]*h*l[7]);
  }
  {
    olim6_mp0 o {{2, 2, 2}, h, S};
    o.add_boundary_node(vec3<double> {x, y, z}, s);
    ASSERT_DOUBLE_EQ(o.get_value({0, 0, 0}), (s + S[0])*h*l[0]/2);
    ASSERT_DOUBLE_EQ(o.get_value({1, 0, 0}), (s + S[1])*h*l[1]/2);
    ASSERT_DOUBLE_EQ(o.get_value({0, 1, 0}), (s + S[2])*h*l[2]/2);
    ASSERT_DOUBLE_EQ(o.get_value({1, 1, 0}), (s + S[3])*h*l[3]/2);
    ASSERT_DOUBLE_EQ(o.get_value({0, 0, 1}), (s + S[4])*h*l[4]/2);
    ASSERT_DOUBLE_EQ(o.get_value({1, 0, 1}), (s + S[5])*h*l[5]/2);
    ASSERT_DOUBLE_EQ(o.get_value({0, 1, 1}), (s + S[6])*h*l[6]/2);
    ASSERT_DOUBLE_EQ(o.get_value({1, 1, 1}), (s + S[7])*h*l[7]/2);
  }
}
