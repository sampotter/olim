#ifndef __NODE_3D_HPP__
#define __NODE_3D_HPP__

#include "abstract_node.hpp"

struct node_3d: public abstract_node {
  static node_3d make_boundary_node(int i, int j, int k, double value = 0.0);
  inline int get_i() const { return _i; }
  inline void set_i(int i) { _i = i; }
  inline int get_j() const { return _j; }
  inline void set_j(int j) { _j = j; }
  inline int get_k() const { return _k; }
  inline void set_k(int k) { _k = k; }
private:
  int _i;
  int _j;
  int _k;
};

#endif // __NODE_3D_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
