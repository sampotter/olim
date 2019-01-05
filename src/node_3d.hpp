#ifndef __NODE_3D_HPP__
#define __NODE_3D_HPP__

#include <src/config.hpp>

#include "abstract_node.hpp"

using index3_t = std::tuple<int, int, int>;

inline
std::string
to_string(index3_t const & index)
{
  std::ostringstream os;
  os << "(" << std::get<0>(index) << ", "
     << std::get<1>(index) << ", " << std::get<2>(index) << ")";
  return os.str();
}

inline
std::ostream &
operator<<(std::ostream & o, index3_t const & index)
{
  return o << "index3_t " << to_string(index);
}

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
  friend std::string to_string(node_3d const & n);
};

inline
std::string
to_string(node_3d const & n)
{
  std::ostringstream os;
  os << to_string(static_cast<abstract_node>(n))
     << ", (" << n.get_i() << ", " << n.get_j() << ", " << n.get_k() << ")";
  return os.str();
}

inline
std::ostream &
operator<<(std::ostream & o, node_3d const & n)
{
  return o << "node_3d {" << to_string(n) << "}";
}

#endif // __NODE_3D_HPP__
