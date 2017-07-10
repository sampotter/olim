#include "heap.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>

#define VALUE(pos) (_data[pos]->get_value())

static int get_left(int pos) {
  return 2*pos + 1;
}

static int get_right(int pos) {
  return 2*pos + 2;
}

static int get_parent(int pos) {
  return (pos - 1)/2;
}

heap::heap(size_t capacity):
  _data {new node *[capacity]},
  _capacity {capacity} {}

heap::~heap() { delete[] _data; }

node *& heap::front() {
  assert(_size > 0);
  return _data[0];
}

node * const & heap::front() const {
  assert(_size > 0);
  return _data[0];
}

bool heap::empty() const {
  return _size == 0;
}

node ** heap::data() const {
  return _data;
}

size_t heap::size() const {
  return _size;
}

void heap::pop_front() {
  swap(0, _size - 1);
  --_size;
  sink(0);
  assert(has_heap_prop());
}

void heap::insert(node * n) {
  if (_size == _capacity) grow();
  n->set_heap_pos(_size);
  assert(_size < _capacity);
  _data[_size++] = n;
  swim(n);
  assert(has_heap_prop());
}

void heap::swim(node * n) {
  assert(_data[n->get_heap_pos()] == n);
  swim(n->get_heap_pos());
}

void heap::swim(int pos) {
  assert(pos < static_cast<int>(_size));
  int parent = get_parent(pos);
  while (pos > 0 && VALUE(parent) > VALUE(pos)) {
    swap(parent, pos);
    pos = parent;
    parent = get_parent(pos);
  }
  assert(has_heap_prop());
}

void heap::sink(int pos) {
  int ch = get_left(pos), n = static_cast<size_t>(size());
  while (ch < n) {
    if (ch + 1 < n && VALUE(ch) > VALUE(ch + 1)) {
      ++ch;
    }
    if (VALUE(pos) > VALUE(ch)) {
      swap(pos, ch);
    }
    pos = ch;
    ch = get_left(pos);
  }
  assert(has_heap_prop());
}

void heap::print() const {
  std::cout << std::endl << "HEAP:" << std::endl;;
  int level = 0;
  int i0 = 0, i1 = 1ul << level;
  while (i0 < static_cast<int>(_size)) {
    std::cout << level << ":";
    for (int i = i0; i < std::min(static_cast<int>(_size), i1); ++i) {
      std::cout << " " << VALUE(i);
    }
    std::cout << std::endl;
    i0 = i1;
    i1 += 1ul << ++level;
  }
  std::cout << std::endl;
}

void heap::grow() {
  _capacity *= 2;
  node ** tmp = new node *[_capacity];
  memcpy(tmp, _data, _size*sizeof(node *));
  delete[] _data;
  _data = tmp;
  assert(has_heap_prop());
}

void heap::swap(int pos1, int pos2) {
	std::swap(_data[pos1], _data[pos2]);
	_data[pos1]->set_heap_pos(pos1);
	_data[pos2]->set_heap_pos(pos2);
}

bool heap::has_heap_prop() const {
  return has_heap_prop_impl(0);
}

bool heap::has_heap_prop_impl(int pos) const {
  double val = VALUE(pos);
  int l = get_left(pos), r = get_right(pos);
  return
    (static_cast<size_t>(l) < size() ? val <= VALUE(l) &&
     has_heap_prop_impl(l) : true) &&
    (static_cast<size_t>(r) < size() ? val <= VALUE(r) &&
     has_heap_prop_impl(r) : true);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
