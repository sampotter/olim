#ifndef __SMART_NODE_HPP__
#define __SMART_NODE_HPP__

#include "abstract_node.hpp"
#include "abstract_smart_node.hpp"

struct smart_node: public abstract_node, public abstract_smart_node {
  struct parent_nodes {
    abstract_node * first {nullptr};
    abstract_node * second {nullptr};
  };
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
