#ifndef __GRAPH_MARCHER_HPP__

#include <vector>

#include "abstract_marcher.hpp"
#include "digraph.hpp"

template <class Node, class Neighbors>
struct graph_marcher: public abstract_marcher {
  using digraph_type = digraph<Node *, Neighbors>;
  
  // TODO: add emplace_back add_node with perfect forwarding
  void add_node(Node const & node);
  void add_neighbor(Node * node, Node * neighbor);

protected:
  Node & get_node(int i) { return _nodes[i]; }
  Node const & get_node(int i) const { return _nodes[i]; }
  void stage(Node * node);
  void update(Node * node);
  
  virtual double speed(Node * node) const = 0;
  virtual Neighbors get_neighbors(abstract_node * node) const;
  virtual void get_valid_neighbors(abstract_node * node, abstract_node ** nb) = 0;
  virtual void update_impl(Node * node, double & T) = 0;
  
private:
  virtual void stage_neighbors_impl(abstract_node * node) override final;
  
  std::vector<Node> _nodes;
  digraph_type _digraph;
};

#include "graph_marcher.impl.hpp"

#endif // __GRAPH_MARCHER_HPP__
