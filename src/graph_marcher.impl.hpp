#ifndef __GRAPH_MARCHER_IMPL_HPP__
#define __GRAPH_MARCHER_IMPL_HPP__

template <class Node, class Neighbors>
void graph_marcher<Node, Neighbors>::add_node(Node const & node) {
  _nodes.push_back(node);
}

template <class Node, class Neighbors>
void graph_marcher<Node, Neighbors>::add_neighbor(Node * node, Node * neighbor) {
  _digraph.add_arc(node, neighbor);
}

template <class Node, class Neighbors>
void graph_marcher<Node, Neighbors>::stage(Node * node) {
  if (node->is_far()) {
    node->set_trial();
    insert_into_heap(node);
  }
}

template <class Node, class Neighbors>
void graph_marcher<Node, Neighbors>::update(Node * node) {
  double T {std::numeric_limits<double>::infinity()};
  update_impl(node, T);
  assert(node->is_trial());
  if (T <= node->get_value()) {
    node->set_value(T);
    adjust_heap_entry(node);
  }
}

template <class Node, class Neighbors>
Neighbors graph_marcher<Node, Neighbors>::get_neighbors(abstract_node * node) const {
  return _digraph.get_neighbors(static_cast<Node *>(node));
}

template <class Node, class Neighbors>
void graph_marcher<Node, Neighbors>::stage_neighbors_impl(abstract_node * node) {
  for (auto neighbor: get_neighbors(node)) {
    stage(neighbor);
  }
  for (auto const & neighbor: get_neighbors(node)) {
    if (!neighbor->is_valid()) {
      update(neighbor);
    }
  }
}

#endif // __GRAPH_MARCHER_IMPL_HPP__
