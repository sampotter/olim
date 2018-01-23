#include <gtest/gtest.h>

#include "heap.hpp"
#include "node.hpp"

TEST (heap, empty_works) {
  heap<node> h {16};
  ASSERT_TRUE(h.empty());
}

TEST (heap, size_works) {
  heap<node> h {16};
  node n[] = {{0, 0, 1}, {0, 0, 2}, {0, 0, 3}};
  ASSERT_TRUE(h.size() == 0);
  h.insert(&n[1]);
  h.insert(&n[0]);
  h.insert(&n[2]);
  ASSERT_TRUE(h.size() == 3);
  h.pop_front();
  ASSERT_TRUE(h.size() == 2);
  h.pop_front();
  h.pop_front();
  ASSERT_TRUE(h.size() == 0);
}

TEST (heap, insert_works) {
  heap<node> h {16};
  node n[] = {{0, 0, 1}, {0, 0, 2}, {0, 0, 3}};
  h.insert(&n[2]);
  h.insert(&n[1]);
  h.insert(&n[0]);
  node * front = h.front();
  ASSERT_TRUE(front == &n[0]);
  ASSERT_TRUE(front->get_value() == 1);
  h.pop_front();
  front = h.front();
  ASSERT_TRUE(front == &n[1]);
  ASSERT_TRUE(front->get_value() == 2);
  h.pop_front();
  front = h.front();
  ASSERT_TRUE(front == &n[2]);
  ASSERT_TRUE(front->get_value() == 3);
  h.pop_front();
}

TEST (heap, swim_works) {
  heap<node> h {16};
  node n[] = {{0, 0, 1}, {0, 0, 2}, {0, 0, 3}};
  h.insert(&n[0]);
  h.insert(&n[1]);
  h.insert(&n[2]);
  n[2].set_value(0);
  h.swim(&n[2]);
  node * front = h.front();
  ASSERT_TRUE(front == &n[2]);
  ASSERT_TRUE(front->get_value() == 0);
}
