#include <gtest/gtest.h>

#include "vec.hpp"

TEST (vec, row_major_indexing_works) {
  {
    vec2<int> inds = {2, 3}, dims = {5, 4};
    int lin = to_linear_index<ordering::ROW_MAJOR>(inds, dims);
    ASSERT_EQ(lin, 11);
    ASSERT_EQ(inds, to_vector_index<ordering::ROW_MAJOR>(lin, dims));
  }
  {
    vec2<int> inds = {5, 4}, dims = {7, 5};
    int lin = to_linear_index<ordering::ROW_MAJOR>(inds, dims);
    ASSERT_EQ(lin, 29);
    ASSERT_EQ(inds, to_vector_index<ordering::ROW_MAJOR>(lin, dims));
  }
  {
    vec3<int> inds = {2, 1, 1}, dims = {4, 3, 4};
    int lin = to_linear_index<ordering::ROW_MAJOR>(inds, dims);
    ASSERT_EQ(lin, 29);
    ASSERT_EQ(inds, to_vector_index<ordering::ROW_MAJOR>(lin, dims));
  }
}

TEST (vec, column_major_indexing_works) {
  {
    vec2<int> inds = {1, 2}, dims = {4, 3};
    int lin = to_linear_index<ordering::COLUMN_MAJOR>(inds, dims);
    ASSERT_EQ(lin, 9);
    ASSERT_EQ(inds, to_vector_index<ordering::COLUMN_MAJOR>(lin, dims));
  }
  {
    vec3<int> inds = {0, 1, 2}, dims = {5, 4, 3};
    int lin = to_linear_index<ordering::COLUMN_MAJOR>(inds, dims);
    ASSERT_EQ(lin, 45);
    ASSERT_EQ(inds, to_vector_index<ordering::COLUMN_MAJOR>(lin, dims));
  }
}
