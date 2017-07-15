#ifndef __NODE_HPP__
#define __NODE_HPP__

#include "abstract_node.hpp"

struct node: public abstract_node {
  static node make_boundary_node(int i, int j, double value = 0.0);
  inline int get_i() const { return _i; }
  inline void set_i(int i) { _i = i; }
  inline int get_j() const { return _j; }
  inline void set_j(int j) { _j = j; }
private:
  int _i;
  int _j;
};

#endif // __NODE_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
