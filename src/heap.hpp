#pragma once

#include <src/config.hpp>

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

#include <algorithm>
#include <iostream>

#include "common.hpp"

inline size_t initial_heap_capacity(int size) {
  return static_cast<size_t>(fmax(8.0, log(size)));
}

template <class elt, class proxy>
struct heap
{
  using value_t = typename proxy::value_t;
  
  heap(proxy && p, size_t capacity);
  ~heap();

  elt & front();
  elt const & front() const;
  bool empty() const;
  elt * data() const;
  size_t size() const;
  void pop_front();
  void insert(elt const & e);
  void update(elt const & e);
  void print() const;

  void grow();
  void swim(int pos);
  void sink(int pos);
  bool has_heap_prop() const;
  bool has_heap_prop_impl(int pos) const;
  void swap(int pos1, int pos2);

  inline value_t value(int pos) const {
    return _proxy.get_value(_data[pos]);
  }
	
  inline int left(int pos) const {
    return 2*pos + 1;
  }

  inline int right(int pos) const {
    return 2*pos + 2;
  }

  inline int parent(int pos) const {
    return (pos - 1)/2;
  }

  proxy _proxy;
  elt * _data {nullptr};
  size_t _size {0};
  size_t _capacity {0};
};

template <class elt, class proxy>
heap<elt, proxy>::heap(proxy && p, size_t capacity):
  _proxy {std::move(p)},
  _data {new elt[capacity]},
  _capacity {capacity}
{}

template <class elt, class proxy>
heap<elt, proxy>::~heap() {
  delete[] _data;
}

template <class elt, class proxy>
elt & heap<elt, proxy>::front() {
  assert(_size > 0);
  return _data[0];
}

template <class elt, class proxy>
elt const & heap<elt, proxy>::front() const {
  assert(_size > 0);
  return _data[0];
}

template <class elt, class proxy>
bool heap<elt, proxy>::empty() const {
  return _size == 0;
}

template <class elt, class proxy>
elt * heap<elt, proxy>::data() const {
  return _data;
}

template <class elt, class proxy>
size_t heap<elt, proxy>::size() const {
  return _size;
}

template <class elt, class proxy>
void heap<elt, proxy>::pop_front() {
  swap(0, _size - 1);
  --_size;
  sink(0);
#if CHECK_HEAP_PROP_IN_DEBUG && OLIM_DEBUG && !RELWITHDEBINFO
  assert(has_heap_prop());
#endif
}

template <class elt, class proxy>
void heap<elt, proxy>::insert(elt const & e) {
  if (_size == _capacity) grow();
  // n->set_heap_pos(_size);
  _proxy.set_heap_pos(e, _size);
  assert(_size < _capacity);
  _data[_size++] = e;
  update(e);
#if CHECK_HEAP_PROP_IN_DEBUG && OLIM_DEBUG && !RELWITHDEBINFO
  assert(has_heap_prop());
#endif
}

template <class elt, class proxy>
void heap<elt, proxy>::update(elt const & e) {
#if OLIM_DEBUG && !RELWITHDEBINFO
  auto const heap_pos = _proxy.get_heap_pos(e);
  assert(_data[heap_pos] == e);
#endif
  swim(_proxy.get_heap_pos(e));
}

template <class elt, class proxy>
void heap<elt, proxy>::swim(int pos) {
  assert(pos < static_cast<int>(_size));
  int parent = this->parent(pos);
  while (pos > 0 && this->value(parent) > this->value(pos)) {
    swap(parent, pos);
    pos = parent;
    parent = this->parent(pos);
  }
#if CHECK_HEAP_PROP_IN_DEBUG && OLIM_DEBUG && !RELWITHDEBINFO
  assert(has_heap_prop());
#endif
}

template <class elt, class proxy>
void heap<elt, proxy>::sink(int pos) {
  int ch = this->left(pos), next = ch + 1, n = static_cast<size_t>(size());
  while (ch < n) {
    if (next < n && this->value(ch) > this->value(next)) {
      ch = next;
    }
    if (this->value(pos) > this->value(ch)) {
      swap(pos, ch);
    }
    pos = ch;
    ch = this->left(pos);
    next = ch + 1;
  }
#if CHECK_HEAP_PROP_IN_DEBUG && OLIM_DEBUG && !RELWITHDEBINFO
  assert(has_heap_prop());
#endif
}

template <class elt, class proxy>
void heap<elt, proxy>::print() const {
  std::cout << std::endl << "HEAP:" << std::endl;;
  int level = 0;
  int i0 = 0, i1 = 1ul << level;
  while (i0 < static_cast<int>(_size)) {
    std::cout << level << ":";
    for (int i = i0; i < std::min(static_cast<int>(_size), i1); ++i) {
      std::cout << " " << this->value(i);
    }
    std::cout << std::endl;
    i0 = i1;
    i1 += 1ul << ++level;
  }
  std::cout << std::endl;
}

template <class elt, class proxy>
void heap<elt, proxy>::grow() {
  _capacity *= 2;
  elt * tmp = new elt[_capacity];
  memcpy(tmp, _data, _size*sizeof(elt));
  delete[] _data;
  _data = tmp;
#if CHECK_HEAP_PROP_IN_DEBUG && OLIM_DEBUG && !RELWITHDEBINFO
  assert(has_heap_prop());
#endif
}

template <class elt, class proxy>
void heap<elt, proxy>::swap(int pos1, int pos2) {
  std::swap(_data[pos1], _data[pos2]);
  _proxy.set_heap_pos(_data[pos1], pos1);
  _proxy.set_heap_pos(_data[pos2], pos2);
}

template <class elt, class proxy>
bool heap<elt, proxy>::has_heap_prop() const {
  return has_heap_prop_impl(0);
}

template <class elt, class proxy>
bool heap<elt, proxy>::has_heap_prop_impl(int pos) const {
  double val = this->value(pos);
  int l = this->left(pos), r = this->right(pos);
  return
    (static_cast<size_t>(l) < size() ? val <= this->value(l) &&
     has_heap_prop_impl(l) : true) &&
    (static_cast<size_t>(r) < size() ? val <= this->value(r) &&
     has_heap_prop_impl(r) : true);
}
