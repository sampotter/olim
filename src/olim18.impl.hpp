#ifndef __OLIM18_IMPL_HPP__
#define __OLIM18_IMPL_HPP__

// neighbor order: U, UN, UE, US, UW, N, NE, E, SE, S, SW, W, NW, D,
// DN, DE, DS, DW

template <class Node, class Updates>
int olim18<Node, Updates>::di[] = {
  0, 1, 0, -1, 0, 1, 1, 0, -1, -1, -1, 0, 1, 0, 1, 0, -1, 0
};

template <class Node, class Updates>
int olim18<Node, Updates>::dj[] = {
  0, 0, 1, 0, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 0, 1, 0, -1
};

template <class Node, class Updates>
int olim18<Node, Updates>::dk[] = {
  1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1
};

template <class Node, class Updates>
void olim18<Node, Updates>::get_valid_neighbors(int i, int j, int k,
                                                abstract_node ** nb) {
  int a, b, c;
  for (int l = 0; l < 18; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->is_valid(a, b, c)) {
      nb[l] = &this->operator()(a, b, c);
    }
  }
}

template <class Node, class Updates>
void olim18<Node, Updates>::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<Node *>(n)->get_i();
  int j = static_cast<Node *>(n)->get_j();
  int k = static_cast<Node *>(n)->get_k();

  for (int l = 0; l < 18; ++l) {
    this->stage(i + di[l], j + dj[l], k + dk[l]);
  }

  int a, b, c;
  for (int l = 0; l < 18; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      this->update(a, b, c);
    }
  }
}

#endif // __OLIM18_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
