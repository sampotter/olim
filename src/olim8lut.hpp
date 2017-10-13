#ifndef __OLIM8LUT_HPP__
#define __OLIM8LUT_HPP__

#include "moore_marcher.hpp"
#include "node.hpp"

template <class update_rules>
struct olim8lut: public moore_marcher<node>, public update_rules {
  using moore_marcher::moore_marcher;
private:
  virtual void update_impl(int i, int j, double & T);
};

#include "olim8lut.impl.hpp"

#endif // __OLIM8LUT_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
