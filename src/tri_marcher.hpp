#ifndef __TRI_MARCHER_HPP__
#define __TRI_MARCHER_HPP__

// TODO: for now, just to get the concept working and to start moving
// forward with the smart marchers, we'll implement this for R2 (i.e.,
// by specializing the node type to graph_node). Later, we'll make
// this more general, and allow for an arbitrary Node type. This will
// take a little more effort...

#include <vector>

#include "graph_node.hpp"

struct tri_marcher;

struct tri_neighbors {
  using iterator = typename std::vector<graph_node *>::iterator;
  void add_neighbor(graph_node * parent, graph_node * neighbor);
  iterator begin();
  iterator end();
private:
  std::vector<graph_node *> _neighbors;
};

struct tri_marcher: public graph_marcher<graph_node, tri_neighbors> {
};

#include "tri_marcher.impl.hpp"

#endif // __TRI_MARCHER_HPP__
