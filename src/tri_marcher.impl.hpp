#ifndef __TRI_MARCHER_IMPL_HPP__
#define __TRI_MARCHER_IMPL_HPP__

#include <cmath>

void tri_neighbors::add_neighbor(graph_node * parent, graph_node * neighbor) {
  double const x0 = parent->get_x();
  double const y0 = parent->get_y();
  double const dx = neighbor->get_x() - x0;
  double const dy = neighbor->get_y() - y0;
  double const theta = atan2(dy, dx);

  // Insert the new neighbor in decreasing order of angle (clockwise
  // order, in keeping with previous implementations)
  //
  // TODO: for now, this just uses insertion sort, since the
  // neighborhood size is likely to be very small. Later, it may be
  // worth conditionally using a heavier duty algorithm here (i.e.,
  // when the number of neighbors exceeds a threshold).
  auto it = _neighbors.begin();
  while (atan2(it->get_y() - y0, it->get_x() - x0) > theta) {
    ++it;
  }
  _neighbors.insert(it, neighbor);
}

tri_neighbors::iterator begin() {
  return _neighbors.begin();
}

tri_neighbors::iterator end() {
  return _neighbors.end();
}

#endif // __TRI_MARCHER_IMPL_HPP__
