#ifndef __NEUMANN_MARCHER_HPP__
#define __NEUMANN_MARCHER_HPP__

#include "marcher.hpp"

template <class Node>
struct neumann_marcher: public marcher<Node> {
  using marcher<Node>::marcher;
EIKONAL_PROTECTED:
  virtual void get_valid_neighbors(int i, int j, abstract_node ** nb);
EIKONAL_PRIVATE:
  virtual void visit_neighbors_impl(abstract_node * n);
};

#include "neumann_marcher.impl.hpp"

#endif // __NEUMANN_MARCHER_HPP__
