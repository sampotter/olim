#ifndef __DIGRAPH_HPP__
#define __DIGRAPH_HPP__

#include <unordered_map>
#include <vector>

template <class V>
struct digraph {
  using neighbors_type = std::vector<V>;
  
  void add_vertex(V const & v);
  void add_arc(V const & u, V const & v);
  neighbors_type const & get_neighbors(V const & v) const;
private:
  std::unordered_map<V, neighbors_type> _adj_list;
};

#include "digraph.impl.hpp"

#endif // __DIGRAPH_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
