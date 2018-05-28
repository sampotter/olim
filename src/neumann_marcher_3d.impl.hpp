#ifndef __NEUMANN_MARCHER_3D_IMPL_HPP__
#define __NEUMANN_MARCHER_3D_IMPL_HPP__

#include <src/config.hpp>

#if PRINT_UPDATES
#    include <cstdio>
#endif

// neighbor order: N, E, U, S, W, D (up -> clockwise from north -> down)

template <class Node>
int neumann_marcher_3d<Node>::di[] = {1, 0, 0, -1, 0, 0};

template <class Node>
int neumann_marcher_3d<Node>::dj[] = {0, 1, 0, 0, -1, 0};

template <class Node>
int neumann_marcher_3d<Node>::dk[] = {0, 0, 1, 0, 0, -1};

template <class Node>
void neumann_marcher_3d<Node>::stage_neighbors_impl(abstract_node * n) {
  int i = static_cast<Node *>(n)->get_i();
  int j = static_cast<Node *>(n)->get_j();
  int k = static_cast<Node *>(n)->get_k();
#if PRINT_UPDATES
  printf("neumann_marcher_3d::stage_neighbors_impl(i = %d, j = %d, k = %d)\n",
         i, j, k);
#endif

  for (int l = 0; l < 6; ++l) {
    this->stage(i + di[l], j + dj[l], k + dk[l]);
  }

  int a, b, c;
  for (int l = 0; l < 6; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      this->update(a, b, c, (l + 3) % 6);
    }
  }
}

template <class Node>
void neumann_marcher_3d<Node>::get_valid_neighbors(int i, int j, int k,
                                                   abstract_node ** nb) {
  int a, b, c;
  for (int l = 0; l < 6; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->in_bounds(a, b, c) && this->is_valid(a, b, c)) {
      nb[l] = &this->operator()(a, b, c);
    }
  }
}

#endif // __NEUMANN_MARCHER_3D_IMPL_HPP__
