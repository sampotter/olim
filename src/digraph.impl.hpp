#ifndef __DIGRAPH_IMPL_HPP__
#define __DIGRAPH_IMPL_HPP__

template <class V>
void digraph<V>::add_vertex(V const & v) {
  _adj_list[v] = std::vector<V>();
}

template <class V>
void digraph<V>::add_arc(V const & u, V const & v) {
  if (_adj_list.find(u) == _adj_list.end()) {
    add_vertex(u);
  }
  // TODO: probably want to implement this so that we insert the
  // nodes sorted into an order which makes sense for our updates
  _adj_list[u].push_back(v);
}

template <class V>
typename digraph<V>::neighbors_type const &
digraph<V>::get_neighbors(V const & v) const {
  return _adj_list.at(v);
}

#endif // __DIGRAPH_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
