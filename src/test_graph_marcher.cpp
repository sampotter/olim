#include "test_graph_marcher.hpp"

#include <cmath>
#include <limits>

test_graph_marcher::test_graph_marcher(int height, int width, double h,
                                       speed_func S, double x0, double y0):
  _S {S}, _height {height}, _width {width}, _h {h}
{
  for (int i = 0; i < _height; ++i) {
    for (int j = 0; j < _width; ++j) {
      add_node({
        h*j - x0,
        h*i - y0,
        std::numeric_limits<double>::infinity(),
        state::far
      });
    }
  }

  auto const in_bounds = [&] (int i, int j) {
    return 0 <= i && i < _height && 0 <= j && j < _width;
  };

  auto const add_neumann_neighborhood = [&] (int i, int j) {
    auto n = &get_node(_width*i + j);
    int a, b;

    a = i - 1, b = j;
    if (in_bounds(a, b)) add_neighbor(n, &get_node(_width*a + b));

    a = i, b = j + 1;
    if (in_bounds(a, b)) add_neighbor(n, &get_node(_width*a + b));

    a = i + 1, b = j;
    if (in_bounds(a, b)) add_neighbor(n, &get_node(_width*a + b));

    a = i, b = j - 1;
    if (in_bounds(a, b)) add_neighbor(n, &get_node(_width*a + b));
  };

  for (int i = 0; i < _height; ++i) {
    for (int j = 0; j < _width; ++j) {
      add_neumann_neighborhood(i, j);
    }
  }
}

void test_graph_marcher::add_boundary_node(int i, int j, double value) {
  auto n = (abstract_node *) &get_node(_width*i + j);
  n->set_valid();
  n->set_value(value);
  stage_neighbors(n);
}

double test_graph_marcher::get_value(int i, int j) const {
  return get_node(_width*i + j).get_value();
}

double test_graph_marcher::S(graph_node * node) const {
  return _S(node->get_x(), node->get_y());
}

void test_graph_marcher::get_valid_neighbors(abstract_node * node,
                                             abstract_node ** nb) {
  int i = 0;
  for (graph_node * neighbor: get_neighbors(node)) {
    nb[i++] = neighbor->is_valid() ? (abstract_node *) neighbor : nullptr;
  }
}

void test_graph_marcher::update_impl(graph_node * node, double & T) {
  abstract_node * nb[4] = {nullptr, nullptr, nullptr, nullptr};
  get_valid_neighbors(node, nb);
  double sh = get_h()*S(node);

  double T1 = std::min(
    nb[0] ? nb[0]->get_value() : std::numeric_limits<double>::infinity(),
    nb[2] ? nb[2]->get_value() : std::numeric_limits<double>::infinity());

  double T2 = std::min(
    nb[1] ? nb[1]->get_value() : std::numeric_limits<double>::infinity(),
    nb[3] ? nb[3]->get_value() : std::numeric_limits<double>::infinity());

  bool T1_inf = std::isinf(T1), T2_inf = std::isinf(T2);

  if (!T1_inf && !T2_inf) {
    double diff = T1 - T2, disc = 2*sh*sh - diff*diff;
    T = disc > 0 ? std::min(T, (T1 + T2 + std::sqrt(disc))/2) : T;
  } else if (std::isinf(T1)) {
    T = std::min(T, T2 + sh);
  } else if (std::isinf(T2)) {
    T = std::min(T, T1 + sh);
  } else {
    assert(false);
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
