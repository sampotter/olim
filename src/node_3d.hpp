#ifndef __NODE_3D_HPP__
#define __NODE_3D_HPP__

#include <src/config.hpp>

#include "abstract_node.hpp"

using index3_t = std::tuple<int, int, int>;

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
#if TRACK_PARENTS
  auto p = n.get_parents();
  node_3d * p0 = static_cast<node_3d *>(p[0]),
    * p1 = static_cast<node_3d *>(p[1]), * p2 = static_cast<node_3d *>(p[2]);
  os << ", parents: ["
     << "(" << p0->get_i() << ", " << p0->get_j() << ", " << p0->get_k() << "), "
     << "(" << p1->get_i() << ", " << p1->get_j() << ", " << p1->get_k() << "), "
     << "(" << p2->get_i() << ", " << p2->get_j() << ", " << p2->get_k() << ")"
     << "]";
#endif
};

#endif // __NODE_3D_HPP__
