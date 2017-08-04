#ifndef __GRAPH_NODE_HPP__
#define __GRAPH_NODE_HPP__

#include "abstract_node.hpp"

struct graph_node: public abstract_node {
  graph_node(double x, double y, double value = 0, state s = state::valid):
    abstract_node {value, s}, _x {x}, _y {y} {}
  inline double get_x() const { return _x; }
  inline double get_y() const { return _y; }
private:
  double _x, _y;
};

#endif // __GRAPH_NODE_HPP_

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
