#ifndef __BASIC_MARCHER_HPP__
#define __BASIC_MARCHER_HPP__

#include "marcher.hpp"
#include "node.hpp"

template <class base>
struct basic_marcher: public base {
  static constexpr int nneib = 4;

EIKONAL_PRIVATE:
  virtual void update_impl(int i, int j, double & T);
};

#include "basic_marcher.impl.hpp"

#endif // __BASIC_MARCHER_HPP__
