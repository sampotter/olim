#ifndef __DIGRAPH_IMPL_HPP__
#define __DIGRAPH_IMPL_HPP__

template <class Vertex, class Neighbors>
void digraph<Vertex, Neighbors>::add_vertex(Vertex const & v) {
  _adj_list[v] = Neighbors {};
}

template <class Vertex, class Neighbors>
void digraph<Vertex, Neighbors>::add_arc(Vertex const & u, Vertex const & v) {
  if (_adj_list.find(u) == _adj_list.end()) {
    add_vertex(u);
  }
  _adj_list[u].add_neighbor(v);
}

template <class Vertex, class Neighbors>
Neighbors const &
digraph<Vertex, Neighbors>::get_neighbors(Vertex const & v) const {
  return _adj_list.at(v);
}

#endif // __DIGRAPH_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
