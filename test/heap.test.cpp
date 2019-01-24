#include <gtest/gtest.h>

#include "heap.hpp"

struct test_proxy
{
  using value_t = int;
  
  test_proxy() {}

  test_proxy(int size):
    _value {new int[size]},
    _heap_pos {new int[size]},
    _size {size}
  {
    memset(_value, 0x0, size*sizeof(int));
    memset(_heap_pos, 0x0, size*sizeof(int));
  }

  test_proxy(int const * v, int size):
    _value {new int[size]},
    _heap_pos {new int[size]},
    _size {size}
  {
    memcpy(_value, v, size*sizeof(int));
    memset(_heap_pos, 0x0, size*sizeof(int));
  }

  test_proxy(test_proxy && p):
    _value {std::exchange(p._value, nullptr)},
    _heap_pos {std::exchange(p._heap_pos, nullptr)},
    _size {std::exchange(p._size, 0)}
  {}

  test_proxy(test_proxy const &) = delete;
  test_proxy & operator=(test_proxy const &) = delete;
  
  ~test_proxy() {
    delete[] _value;
    delete[] _heap_pos;
  }

  inline int get_heap_pos(int lin) const {
    return _heap_pos[lin];
  }

  inline void set_heap_pos(int lin, int pos) {
    _heap_pos[lin] = pos;
  }

  inline value_t get_value(int lin) const {
    return _value[lin];
  }

  value_t * _value {nullptr};
  int * _heap_pos {nullptr};
  int _size {0};
};

TEST (heap, empty_works) {
  heap<int, test_proxy> h {{}, 16};
  ASSERT_TRUE(h.empty());
}

TEST (heap, size_works) {
  int const values[] = {1, 2, 3};
  heap<int, test_proxy> h {{values, 3}, 16};
  ASSERT_TRUE(h.size() == 0);
  h.insert(1);
  h.insert(0);
  h.insert(2);
  ASSERT_TRUE(h.size() == 3);
  h.pop_front();
  ASSERT_TRUE(h.size() == 2);
  h.pop_front();
  h.pop_front();
  ASSERT_TRUE(h.size() == 0);
}

TEST (heap, insert_works) {
  int const values[] = {1, 2, 3};
  heap<int, test_proxy> h {{values, 3}, 16};
  h.insert(2);
  h.insert(1);
  h.insert(0);
  int front = h.front();
  ASSERT_TRUE(front == 0);
  ASSERT_TRUE(values[front] == 1);
  h.pop_front();
  front = h.front();
  ASSERT_TRUE(front == 1);
  ASSERT_TRUE(values[front] == 2);
  h.pop_front();
  front = h.front();
  ASSERT_TRUE(front == 2);
  ASSERT_TRUE(values[front] == 3);
  h.pop_front();
}

TEST (heap, update_works) {
  int const values[] = {1, 2, 3};
  heap<int, test_proxy> h {{values, 3}, 16};
  h.insert(0);
  h.insert(1);
  h.insert(2);
  h._proxy._value[2] = 0;
  h.update(2);
  int front = h.front();
  ASSERT_TRUE(front == 2);
  ASSERT_TRUE(h._proxy._value[front] == 0);
}
