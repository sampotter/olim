#ifndef __SMART_NEUMANN_MARCHER_HPP__
#define __SMART_NEUMANN_MARCHER_HPP__

#include "smart_marcher.hpp"

struct smart_neumann_marcher: public smart_marcher {
  using smart_marcher::smart_marcher;
protected:
  virtual void get_valid_neighbors(int i, int j, abstract_node ** nb);
  static int di[4];
  static int dj[4];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
};

#endif // __SMART_NEUMANN_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
