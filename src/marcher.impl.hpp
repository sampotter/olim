#ifndef __MARCHER_IMPL_HPP_HPP__
#define __MARCHER_IMPL_HPP_HPP__

#include <cassert>
#include <cmath>

static inline size_t get_initial_heap_size(int width, int height) {
  return static_cast<size_t>(std::max(8.0, std::log(width*height)));
}

template <class Node>
marcher<Node>::marcher(int height, int width, double h,
                       std::function<double(double, double)> s,
                       double x0, double y0):
  abstract_marcher {get_initial_heap_size(width, height)},
  _nodes {new Node[width*height]},
  _s_cache {new double[width*height]},
  _h {h},
  _x0 {x0},
  _y0 {y0},
  _height {height},
  _width {width}
{
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      _s_cache[i*width + j] = s(h*j - x0, h*i - y0);
    }
  }

  init();
}

template <class Node>
marcher<Node>::marcher(int height, int width, double h,
                       std::unique_ptr<double[]> s_cache):
  abstract_marcher {get_initial_heap_size(width, height)},
  _nodes {new Node[width*height]},
  _s_cache {std::move(s_cache)},
  _h {h},
  _height {height},
  _width {width}
{
  init();
}

template <class Node>
void marcher<Node>::add_boundary_node(int i, int j, double value) {
  assert(in_bounds(i, j));
  assert(operator()(i, j).is_far()); // TODO: for now---worried about heap
  stage_neighbors(&(operator()(i, j) = {i, j, value}));
}

template <class Node>
double marcher<Node>::get_value(int i, int j) const {
  assert(in_bounds(i, j));
  return operator()(i, j).get_value();
}

template <class Node>
Node & marcher<Node>::operator()(int i, int j) {
  assert(in_bounds(i, j));
  return _nodes[_width*i + j];
}

template <class Node>
Node const & marcher<Node>::operator()(int i, int j) const {
  assert(in_bounds(i, j));
  return _nodes[_width*i + j];
}

template <class Node>
void marcher<Node>::stage(int i, int j) {
  if (in_bounds(i, j) && operator()(i, j).is_far()) {
    operator()(i, j).set_trial();
    insert_into_heap(&operator()(i, j));
  }
}

template <class Node>
void marcher<Node>::update(int i, int j) {
  assert(in_bounds(i, j));
  double T = std::numeric_limits<double>::infinity();
  update_impl(i, j, T);
  auto n = &operator()(i, j);
  assert(n->is_trial());
  if (T <= n->get_value()) {
    n->set_value(T);
    adjust_heap_entry(n);
  }
}

template <class Node>
bool marcher<Node>::in_bounds(int i, int j) const {
  return (unsigned) i < (unsigned) _height && (unsigned) j < (unsigned) _width;
}

template <class Node>
double marcher<Node>::speed(int i, int j) {
  assert(in_bounds(i, j));
  return _s_cache[_width*i + j];
}

template <class Node>
bool marcher<Node>::is_valid(int i, int j) const {
  return in_bounds(i, j) && operator()(i, j).is_valid();
}

template <class Node>
void marcher<Node>::init() {
  for (int i = 0; i < _height; ++i) {
    for (int j = 0; j < _width; ++j) {
      operator()(i, j).set_i(i);
      operator()(i, j).set_j(j);
    }
  }
}

#endif // __MARCHER_IMPL_HPP_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
