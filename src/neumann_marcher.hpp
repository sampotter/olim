#ifndef __NEUMANN_MARCHER_HPP__
#define __NEUMANN_MARCHER_HPP__

#include "fast_marcher.hpp"

struct neumann_marcher: public fast_marcher {
  using fast_marcher::fast_marcher;
protected:
  virtual void get_valid_neighbors(int i, int j, abstract_node ** nb);
  static int di[4];
  static int dj[4];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
};

#endif // __NEUMANN_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
