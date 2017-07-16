#ifndef __ABSTRACT_SMART_NODE_HPP__
#define __ABSTRACT_SMART_NODE_HPP__

#include "abstract_smart_node.hpp"

struct abstract_smart_node {
  inline double get_lambda() const { return _lambda; }
  inline void set_lambda(double lambda) { _lambda = lambda; }
private:
  double _lambda {-1};
};

#endif // __ABSTRACT_SMART_NODE_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
