#ifndef __MOORE_MARCHER_IMPL_HPP__
#define __MOORE_MARCHER_IMPL_HPP__

#include <src/config.hpp>

#include "offsets.hpp"

#define __di(k) di<2>[k]
#define __dj(k) dj<2>[k]

template <class Node>
void moore_marcher<Node>::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<Node *>(n)->get_i();
  int j = static_cast<Node *>(n)->get_j();

  // TODO: see comment in neumann_marcher::stage_neighbors_impl

  for (int k = 0; k < 8; ++k) this->pre_stage(i + __di(k), j + __dj(k));

  int a, b;
  for (int k = 0; k < 8; ++k) {
    a = i + __di(k), b = j + __dj(k);
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
    a = i + __di(k), b = j + __dj(k);
    if (this->in_bounds(a, b) && this->is_valid(a, b)) {
      nb[k] = &this->operator()(a, b);
    }
  }
}

#undef __di
#undef __dj

#endif // __MOORE_MARCHER_IMPL_HPP__
