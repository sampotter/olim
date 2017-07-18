#ifndef __SMART_MARCHER_HPP__
#define __SMART_MARCHER_HPP__

#include "marcher.hpp"
#include "smart_node.hpp"

struct smart_marcher: public marcher<smart_node> {
  using marcher<smart_node>::marcher;
  
  parent_nodes get_parent_nodes(int i, int j) const {
    return operator()(i, j).get_parent_nodes();
  }
  
  double get_lambda(int i, int j) const {
    return operator()(i, j).get_lambda();
  };
};

#endif // __SMART_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
