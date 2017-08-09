#ifndef __OLIM8FAST_HPP__
#define __OLIM8FAST_HPP__

#include "moore_marcher.hpp"
#include "node.hpp"

template <class update_rules>
struct olim8fast: public moore_marcher<node>, public update_rules {
  using moore_marcher::moore_marcher;
private:
  virtual void update_impl(int i, int j, double & T);
};

#include "olim8fast.impl.hpp"

#endif // __OLIM8FAST_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
