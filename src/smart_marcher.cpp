#include "smart_marcher.hpp"

smart_marcher::smart_marcher(int height, int width, double h, speed_func S,
                             double x0, double y0):
  marcher {height, width, h, S, x0, y0},
  _nodes {new smart_node[width*height]}
{
  init();
}

smart_marcher::smart_marcher(int height, int width, double h, double * S_cache):
  marcher {height, width, h, S_cache},
  _nodes {new smart_node[width*height]}
{
  init();
}

void smart_marcher::add_boundary_node(int i, int j, double value) {
  assert(in_bounds(i, j));
  assert(operator()(i, j).is_far()); // TODO: for now---worried about heap
  operator()(i, j) = smart_node::make_boundary_node(i, j, value);
  stage_neighbors(&operator()(i, j));
}

double smart_marcher::get_value(int i, int j) const {
  assert(in_bounds(i, j));
  return operator()(i, j).get_value();
}

smart_node & fast_marcher::operator()(int i, int j) {
  assert(in_bounds(i, j));
  return _nodes[get_width()*i + j];
}

smart_node const & fast_marcher::operator()(int i, int j) const {
  assert(in_bounds(i, j));
  return _nodes[get_width()*i + j];
}

void smart_marcher::stage_neighbor(int i, int j) {
  if (in_bounds(i, j) && operator()(i, j).is_far()) {
    operator()(i, j).set_trial();
    insert_into_heap(&operator()(i, j));
  }
}

// TODO: this is where we want to update the smart nodes
void smart_marcher::update_node_value(int i, int j) {
  assert(in_bounds(i, j));
  double T = std::numeric_limits<double>::infinity();
  update_node_value_impl(i, j, T);
  smart_node * n = &operator()(i, j);
  assert(n->is_trial());
  if (T <= n->get_value()) {
    n->set_value(T);
    adjust_heap_entry(n);
  }
}

bool smart_marcher::is_valid(int i, int j) const {
  return in_bounds(i, j) && operator()(i, j).is_valid();
}

void smart_marcher::init() {
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
