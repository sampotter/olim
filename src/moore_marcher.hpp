#ifndef __MOORE_MARCHER_HPP__
#define __MOORE_MARCHER_HPP__

#include "marcher.hpp"

template <class Node>
struct moore_marcher: public marcher<Node> {
  using marcher<Node>::marcher;
protected:
  virtual void get_valid_neighbors(int i, int j, abstract_node ** nb);
private:
  virtual void stage_neighbors_impl(abstract_node * n);
};

#include "moore_marcher.impl.hpp"

#endif // __MOORE_MARCHER_HPP__
