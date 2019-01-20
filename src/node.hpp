#ifndef __NODE_HPP__
#define __NODE_HPP__

#include "abstract_node.hpp"

#include <math.h>

struct node: public abstract_node {
  struct fac_center {
    fac_center(double i, double j, double s): i {i}, j {j}, s {s} {}
    double i, j, s;
  };

  node() {}
  node(int i, int j, double value = 0, state s = state::valid):
    abstract_node {value, s}, _i {i}, _j {j} {}

  inline int get_i() const { return _i; }
  inline void set_i(int i) { _i = i; }
  inline int get_j() const { return _j; }
  inline void set_j(int j) { _j = j; }

  inline bool is_factored() const { return _fac_center != nullptr; }
  inline const fac_center * get_fac_center() const { return _fac_center; }
  inline void set_fac_center(fac_center const * fc) { _fac_center = fc; }

  // TODO: these are candidates for a bit more optimization, probably
  // in the constructor of this class...
  inline int get_i_fac() const {
    return _i - _fac_center->i;
  }
  inline int get_j_fac() const {
    return _j - _fac_center->j;
  }

EIKONAL_PRIVATE:
  int _i {-1}, _j {-1};
  fac_center const * _fac_center {nullptr};
};

#endif // __NODE_HPP__
