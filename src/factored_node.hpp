#ifndef __FACTORED_NODE_HPP__
#define __FACTORED_NODE_HPP__

#include "node.hpp"

struct factored_node: public node {
  factored_node() {}
  factored_node(factored_node * src, double dist, int i, int j, double value = 0):
    node {i, j, value} {}
  inline factored_node * src() const { return _src; }
  inline double dist() const { return _dist; }
EIKONAL_PRIVATE:
  factored_node * _src;
};

#endif // __FACTORED_NODE_HPP__
