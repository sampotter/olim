#ifndef __NEUMANN_MARCHER_IMPL_HPP__
#define __NEUMANN_MARCHER_IMPL_HPP__

// neighbor order: N, E, S, W (clockwise from north)

template <class Node>
int neumann_marcher<Node>::di[] = {-1, 0, 1, 0};

template <class Node>
int neumann_marcher<Node>::dj[] = {0, 1, 0, -1};

template <class Node>
void neumann_marcher<Node>::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<Node *>(n)->get_i();
  int j = static_cast<Node *>(n)->get_j();

  for (int k = 0; k < 4; ++k) {
    this->stage_neighbor(i + di[k], j + dj[k]);
  }

  int a, b;
  for (int k = 0; k < 4; ++k) {
    a = i + di[k], b = j + dj[k];
    if (this->in_bounds(a, b) && !this->operator()(a, b).is_valid()) {
      this->update(a, b);
    }
  }
}

template <class Node>
void neumann_marcher<Node>::get_valid_neighbors(int i, int j,
                                                abstract_node ** nb) {
  int a, b;
  for (int k = 0; k < 4; ++k) {
    a = i + di[k], b = j + dj[k];
    if (this->is_valid(a, b)) {
      nb[k] = &this->operator()(a, b);
    }
  }
}

#endif // __NEUMANN_MARCHER_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
