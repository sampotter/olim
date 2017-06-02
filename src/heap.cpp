#include "heap.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <functional>
#include <iostream>

static size_t get_left(size_t pos) {
  return 2*pos + 1;
}

static size_t get_right(size_t pos) {
  return 2*pos + 2;
}

static size_t get_parent(size_t pos) {
  return (pos - 1)/2;
}

heap::heap(size_t capacity):
  _data {new node*[capacity]},
  _capacity {capacity} {}

heap::~heap() { delete[] _data; }

node* & heap::front() {
  assert(_size > 0);
  // assert(has_heap_property());
  return _data[0];
}

node* const & heap::front() const {
  assert(_size > 0);
  // assert(has_heap_property());
  return _data[0];
}

bool heap::empty() const {
  return _size == 0;
}

node** heap::data() const {
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

void heap::insert(node* n) {
  if (_size == _capacity) grow();
  n->set_heap_pos(_size);
  _data[_size++] = n;
  adjust_entry(n);
}

void heap::adjust_entry(node* n) {
  size_t pos = n->get_heap_pos();
  assert(_data[pos] == n);
  assert(pos < _size);
  size_t parent = get_parent(pos);
  while (pos > 0 && get_value(parent) > get_value(pos)) {
    swap(parent, pos);
    pos = parent;
    parent = get_parent(pos);
  }
  // assert(has_heap_property());
}

void heap::print() const {
  std::cout << std::endl << "HEAP:" << std::endl;;
  size_t level = 0;
  size_t i0 = 0, i1 = 1ul << level;
  while (i0 < _size) {
    std::cout << level << ":";
    for (size_t i = i0; i < std::min(_size, i1); ++i) {
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
  node** tmp = new node*[_capacity];
  memcpy(tmp, _data, _size*sizeof(node*));
  delete[] _data;
  _data = tmp;
}

void heap::heapify(size_t pos) {
  std::function<void(size_t)> const rec = [&] (size_t pos) {
    size_t l = get_left(pos), r = get_right(pos);
    size_t smallest;
    if (l < _size && get_value(l) < get_value(pos)) {
      smallest = l;
    } else {
      smallest = pos;
    }
    if (r < _size && get_value(r) < get_value(smallest)) {
      smallest = r;
    }
    if (smallest != pos) {
      swap(pos, smallest);
      rec(smallest);
    }
  };
  rec(pos);
  // assert(has_heap_property());
}

bool heap::has_heap_property() const {
  if (empty()) {
    return true;
  }
  std::function<bool(size_t)> const rec = [&] (size_t pos) {
    bool has_prop = true;
    size_t const l = get_left(pos), r = get_right(pos);
    if (l < _size) {
      has_prop &= get_value(pos) <= get_value(l);
      has_prop &= rec(l);
    }
    if (r < _size) {
      has_prop &= get_value(pos) <= get_value(r);
      has_prop &= rec(r);
    }
    return has_prop;
  };
return rec(0);
}

void heap::swap(size_t pos1, size_t pos2) {
	std::swap(_data[pos1], _data[pos2]);
	_data[pos1]->set_heap_pos(pos1);
	_data[pos2]->set_heap_pos(pos2);
}

double heap::get_value(size_t pos) const {
	return _data[pos]->get_value();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
