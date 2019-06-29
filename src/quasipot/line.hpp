#pramga once

#include "../vec.hpp"

namespace quasipot {

template <int n>
struct line
{
  using vec_t = vec<double, n>;

  double operator(double U0, vec const & b, vec const & x0) const {
    return U0 + b.norm2()*x0.norm() - b*x0;
  }
};

}
