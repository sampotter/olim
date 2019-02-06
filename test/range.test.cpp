#include <gtest/gtest.h>

#include "range.hpp"

#include <cstdio>

TEST (range, begin_is_correct) {
  range<2> r {{2, 3}};
  auto begin = vec2<int>::zero();
  ASSERT_EQ(*r.begin(), begin);
}

TEST (range, end_is_correct) {
  range<2> r {{2, 3}};
  vec2<int> end {0, 3};
  ASSERT_EQ(*r.end(), end);
}

TEST (range, iteration_is_correct) {
  range<2> r {{2, 3}};
  ASSERT_NE(*r++, *r.end());
  ASSERT_NE(*r++, *r.end());
  ASSERT_NE(*r++, *r.end());
  ASSERT_NE(*r++, *r.end());
  ASSERT_NE(*r++, *r.end());
  ASSERT_NE(*r++, *r.end());
  ASSERT_EQ(*r++, *r.end());
}

TEST (range, ordering_in_2d_is_correct) {
  {
    vec2<int> inds[4] = {
      {0, 0},
      {1, 0},
      {0, 1},
      {1, 1}
    };
    {
      range<2> r {{2, 2}};
      ASSERT_EQ(*r++, inds[0]);
      ASSERT_EQ(*r++, inds[1]);
      ASSERT_EQ(*r++, inds[2]);
      ASSERT_EQ(*r++, inds[3]);
    }
    {
      int i = 0;
      for (auto inds_: range<2> {{2, 2}}) {
        ASSERT_EQ(inds_, inds[i++]);
      }
    }
  }
  {
    vec2<int> inds[6] = {
      {0, 0},
      {1, 0},
      {0, 1},
      {1, 1},
      {0, 2},
      {1, 2}
    };
    int i = 0;
    for (auto inds_: range<2> {{2, 3}}) {
      ASSERT_EQ(inds_, inds[i++]);
    }
  }
}

TEST (range, ordering_in_3d_is_correct) {
  vec3<int> inds[24] = {
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {1, 1, 0},
    {0, 2, 0},
    {1, 2, 0},
    {0, 0, 1},
    {1, 0, 1},
    {0, 1, 1},
    {1, 1, 1},
    {0, 2, 1},
    {1, 2, 1},
    {0, 0, 2},
    {1, 0, 2},
    {0, 1, 2},
    {1, 1, 2},
    {0, 2, 2},
    {1, 2, 2},
    {0, 0, 3},
    {1, 0, 3},
    {0, 1, 3},
    {1, 1, 3},
    {0, 2, 3},
    {1, 2, 3}
  };
  int i = 0;
  for (auto inds_: range<3> {{2, 3, 4}}) {
    ASSERT_EQ(inds_, inds[i++]);
  }
}
