#ifndef __NEUMANN_MARCHER_HPP__
#define __NEUMANN_MARCHER_HPP__

#include "marcher.hpp"

template <class Node>
struct neumann_marcher: public marcher<Node> {
  using marcher<Node>::marcher;
protected:
  virtual void get_valid_neighbors(int i, int j, abstract_node ** nb);
  static int di[4];
  static int dj[4];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
};

#include "neumann_marcher.impl.hpp"

#endif // __NEUMANN_MARCHER_HPP__
