#ifndef __HEAP_HPP__
#define __HEAP_HPP__

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <iostream>

template <class Node>
struct heap
{
  heap(size_t capacity);
  ~heap();

  Node *& front();
  Node * const & front() const;
  bool empty() const;
  Node ** data() const;
  size_t size() const;
  void pop_front();
  void insert(Node * n);
  void swim(Node * n);
  void print() const;
private:
  void grow();
  void swim(int pos);
  void sink(int pos);
  bool has_heap_prop() const;
  bool has_heap_prop_impl(int pos) const;
  void swap(int pos1, int pos2);
	
  Node ** _data {nullptr};
  size_t _size {0};
  size_t _capacity {0};
};

#define VALUE(pos) (_data[pos]->get_value())
#define LEFT(pos) (2*(pos) + 1)
#define RIGHT(pos) (2*(pos) + 2)
#define PARENT(pos) (((pos) - 1)/2)

template <class Node>
heap<Node>::heap(size_t capacity):
  _data {new Node *[capacity]},
  _capacity {capacity}
{}

template <class Node>
heap<Node>::~heap() {
  delete[] _data;
}

template <class Node>
Node *& heap<Node>::front() {
  assert(_size > 0);
  return _data[0];
}

template <class Node>
Node * const & heap<Node>::front() const {
  assert(_size > 0);
  return _data[0];
}

template <class Node>
bool heap<Node>::empty() const {
  return _size == 0;
}

template <class Node>
Node ** heap<Node>::data() const {
  return _data;
}

template <class Node>
size_t heap<Node>::size() const {
  return _size;
}

template <class Node>
void heap<Node>::pop_front() {
  swap(0, _size - 1);
  --_size;
  sink(0);
  assert(has_heap_prop());
}

template <class Node>
void heap<Node>::insert(Node * n) {
  if (_size == _capacity) grow();
  n->set_heap_pos(_size);
  assert(_size < _capacity);
  _data[_size++] = n;
  swim(n);
  assert(has_heap_prop());
}

template <class Node>
void heap<Node>::swim(Node * n) {
  assert(_data[n->get_heap_pos()] == n);
  swim(n->get_heap_pos());
}

template <class Node>
void heap<Node>::swim(int pos) {
  assert(pos < static_cast<int>(_size));
  int parent = PARENT(pos);
  while (pos > 0 && VALUE(parent) > VALUE(pos)) {
    swap(parent, pos);
    pos = parent;
    parent = PARENT(pos);
  }
  assert(has_heap_prop());
}

template <class Node>
void heap<Node>::sink(int pos) {
  int ch = LEFT(pos), n = static_cast<size_t>(size());
  while (ch < n) {
    if (ch + 1 < n && VALUE(ch) > VALUE(ch + 1)) {
      ++ch;
    }
    if (VALUE(pos) > VALUE(ch)) {
      swap(pos, ch);
    }
    pos = ch;
    ch = LEFT(pos);
  }
  assert(has_heap_prop());
}

template <class Node>
void heap<Node>::print() const {
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

template <class Node>
void heap<Node>::grow() {
  _capacity *= 2;
  Node ** tmp = new Node *[_capacity];
  memcpy(tmp, _data, _size*sizeof(Node *));
  delete[] _data;
  _data = tmp;
  assert(has_heap_prop());
}

template <class Node>
void heap<Node>::swap(int pos1, int pos2) {
	std::swap(_data[pos1], _data[pos2]);
	_data[pos1]->set_heap_pos(pos1);
	_data[pos2]->set_heap_pos(pos2);
}

template <class Node>
bool heap<Node>::has_heap_prop() const {
  return has_heap_prop_impl(0);
}

template <class Node>
bool heap<Node>::has_heap_prop_impl(int pos) const {
  double val = VALUE(pos);
  int l = LEFT(pos), r = RIGHT(pos);
  return
    (static_cast<size_t>(l) < size() ? val <= VALUE(l) &&
     has_heap_prop_impl(l) : true) &&
    (static_cast<size_t>(r) < size() ? val <= VALUE(r) &&
     has_heap_prop_impl(r) : true);
}

#undef PARENT
#undef RIGHT
#undef LEFT
#undef VALUE

#endif // __HEAP_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
