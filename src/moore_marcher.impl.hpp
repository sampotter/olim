#ifndef __MOORE_MARCHER_IMPL_HPP__
#define __MOORE_MARCHER_IMPL_HPP__

// neighbor order: N, NE, E, SE, S, SW, W, NW (clockwise from north)

template <class Node>
int moore_marcher<Node>::di[] = {-1, -1, 0, 1, 1, 1, 0, -1};

template <class Node>
int moore_marcher<Node>::dj[] = {0, 1, 1, 1, 0, -1, -1, -1};

template <class Node>
void moore_marcher<Node>::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<Node *>(n)->get_i();
  int j = static_cast<Node *>(n)->get_j();

  for (int k = 0; k < 8; ++k) {
    this->stage_neighbor(i + di[k], j + dj[k]);
  }

  int a, b;
  for (int k = 0; k < 8; ++k) {
    a = i + di[k];
    b = j + dj[k];
    if (this->in_bounds(a, b) && !this->operator()(a, b).is_valid()) {
      this->update(a, b);
    }
  }
}

template <class Node>
void moore_marcher<Node>::get_valid_neighbors(int i, int j,
                                              abstract_node ** nb) {
  int a, b;
  for (int k = 0; k < 8; ++k) {
    a = i + di[k];
    b = j + dj[k];
    if (this->is_valid(a, b)) {
      nb[k] = &this->operator()(a, b);
    }
  }
}

#endif // __MOORE_MARCHER_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
