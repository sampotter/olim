#ifndef __NEUMANN_MARCHER_IMPL_HPP__
#define __NEUMANN_MARCHER_IMPL_HPP__

#include "offsets.hpp"

#define __di(k) di<2>[k]
#define __dj(k) dj<2>[k]

template <class Node>
void neumann_marcher<Node>::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<Node *>(n)->get_i();
  int j = static_cast<Node *>(n)->get_j();

  for (int k = 0; k < 4; ++k) {
    this->stage(i + __di(k), j + __dj(k));
  }

  int a, b;
  for (int k = 0; k < 4; ++k) {
    a = i + __di(k), b = j + __dj(k);
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
    a = i + __di(k), b = j + __dj(k);
    if (this->is_valid(a, b)) {
      nb[k] = &this->operator()(a, b);
    }
  }
}

#undef __di
#undef __dj

#endif // __NEUMANN_MARCHER_IMPL_HPP__
