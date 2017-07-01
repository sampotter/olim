#include "heap.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <functional>
#include <iostream>

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
  heapify(0);
}

void heap::insert(node * n) {
  if (_size == _capacity) grow();
  n->set_heap_pos(_size);
  assert(_size < _capacity);
  _data[_size++] = n;
  adjust_entry(n);
}

void heap::adjust_entry(node * n) {
  int pos = n->get_heap_pos();
  assert(_data[pos] == n);
  assert(pos < static_cast<int>(_size));
  int parent = get_parent(pos);
  while (pos > 0 && get_value(parent) > get_value(pos)) {
    swap(parent, pos);
    pos = parent;
    parent = get_parent(pos);
  }
}

void heap::print() const {
  std::cout << std::endl << "HEAP:" << std::endl;;
  int level = 0;
  int i0 = 0, i1 = 1ul << level;
  while (i0 < static_cast<int>(_size)) {
    std::cout << level << ":";
    for (int i = i0; i < std::min(static_cast<int>(_size), i1); ++i) {
      std::cout << " " << get_value(i);
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
}

void heap::heapify(int pos) {
  std::function<void(int)> const rec = [&] (int pos) {
    int l = get_left(pos), r = get_right(pos);
    int smallest;
    if (l < static_cast<int>(_size) && get_value(l) < get_value(pos)) {
      smallest = l;
    } else {
      smallest = pos;
    }
    if (r < static_cast<int>(_size) && get_value(r) < get_value(smallest)) {
      smallest = r;
    }
    if (smallest != pos) {
      swap(pos, smallest);
      rec(smallest);
    }
  };
  rec(pos);
}

void heap::swap(int pos1, int pos2) {
	std::swap(_data[pos1], _data[pos2]);
	_data[pos1]->set_heap_pos(pos1);
	_data[pos2]->set_heap_pos(pos2);
}

double heap::get_value(int pos) const {
	return _data[pos]->get_value();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
