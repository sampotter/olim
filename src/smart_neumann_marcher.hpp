#ifndef __SMART_NEUMANN_MARCHER_HPP__
#define __SMART_NEUMANN_MARCHER_HPP__

#include "neumann_marcher.hpp"
#include "smart_node.hpp"

struct smart_neumann_marcher: public neumann_marcher<smart_node> {
  using neumann_marcher::neumann_marcher;

  double get_lambda(int i, int j) const {
    return operator()(i, j).get_lambda();
  }

  parent_nodes const & get_parent_nodes(int i, int j) const {
    return operator()(i, j).get_parent_nodes();
  }
};

#endif // __SMART_NEUMANN_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
