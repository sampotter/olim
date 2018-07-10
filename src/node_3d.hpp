#ifndef __NODE_3D_HPP__
#define __NODE_3D_HPP__

#include "abstract_node.hpp"

struct node_3d: public abstract_node {
  using abstract_node::abstract_node;
  node_3d() {}
  node_3d(int i, int j, int k, double value = 0):
    abstract_node {value, state::valid}, _i {i}, _j {j}, _k {k} {}
  inline int get_i() const { return _i; }
  inline void set_i(int i) { _i = i; }
  inline int get_j() const { return _j; }
  inline void set_j(int j) { _j = j; }
  inline int get_k() const { return _k; }
  inline void set_k(int k) { _k = k; }
EIKONAL_PRIVATE:
  int _i {-1}, _j {-1}, _k {-1};
};

#endif // __NODE_3D_HPP__
