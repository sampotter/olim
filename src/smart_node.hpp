#ifndef __SMART_NODE_HPP__
#define __SMART_NODE_HPP__

#include "abstract_smart_node.hpp"
#include "node.hpp"

struct parent_nodes {
  abstract_node * first {nullptr};
  abstract_node * second {nullptr};
};

struct smart_node: public node, public abstract_smart_node {
  using node::node;
  parent_nodes get_parent_nodes() const { return _parent_nodes; }
  void set_parent_nodes(parent_nodes p) { _parent_nodes = p; }
private:
  parent_nodes _parent_nodes;
};

#endif // __SMART_NODE_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
