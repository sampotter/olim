#include "fast_marcher.hpp"

fast_marcher::fast_marcher(int height, int width, double h, speed_func S,
                           double x0, double y0):
  marcher {height, width, h, S, x0, y0},
  _nodes {new node[width*height]}
{
  init();
}

fast_marcher::fast_marcher(int height, int width, double h, double * S_cache):
  marcher {height, width, h, S_cache},
  _nodes {new node[width*height]}
{
  init();
}

void fast_marcher::add_boundary_node(int i, int j, double value) {
  assert(in_bounds(i, j));
  assert(operator()(i, j).is_far()); // TODO: for now---worried about heap
  operator()(i, j) = node::make_boundary_node(i, j, value);
  stage_neighbors(&operator()(i, j));
}

double fast_marcher::get_value(int i, int j) const {
  assert(in_bounds(i, j));
  return operator()(i, j).get_value();
}

node & fast_marcher::operator()(int i, int j) {
  assert(in_bounds(i, j));
  return _nodes[get_width()*i + j];
}

node const & fast_marcher::operator()(int i, int j) const {
  assert(in_bounds(i, j));
  return _nodes[get_width()*i + j];
}

void fast_marcher::stage_neighbor(int i, int j) {
  if (in_bounds(i, j) && operator()(i, j).is_far()) {
    operator()(i, j).set_trial();
    insert_into_heap(&operator()(i, j));
  }
}

void fast_marcher::update_node_value(int i, int j) {
  assert(in_bounds(i, j));
  double T = std::numeric_limits<double>::infinity();
  update_node_value_impl(i, j, T);
  node * n = &operator()(i, j);
  assert(n->is_trial());
  if (T <= n->get_value()) {
    n->set_value(T);
    adjust_heap_entry(n);
  }
}

bool fast_marcher::is_valid(int i, int j) const {
  return in_bounds(i, j) && operator()(i, j).is_valid();
}

void fast_marcher::init() {
  for (int i = 0; i < get_height(); ++i) {
    for (int j = 0; j < get_width(); ++j) {
      operator()(i, j).set_i(i);
      operator()(i, j).set_j(j);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
