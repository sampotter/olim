#ifndef __DIGRAPH_HPP__
#define __DIGRAPH_HPP__

#include <unordered_map>

template <class Vertex, class Neighbors>
struct digraph {
  void add_vertex(Vertex const & v);
  void add_arc(Vertex const & u, Vertex const & v);
  Neighbors const & get_neighbors(Vertex const & v) const;
private:
  std::unordered_map<Vertex, Neighbors> _adj_list;
};

#include "digraph.impl.hpp"

#endif // __DIGRAPH_HPP__
