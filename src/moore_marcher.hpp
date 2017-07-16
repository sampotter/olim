#ifndef __MOORE_MARCHER_HPP__
#define __MOORE_MARCHER_HPP__

#include "marcher.hpp"

struct moore_marcher: public marcher
{
  using marcher::marcher;
protected:
  virtual void get_valid_neighbors(int i, int j, node ** nb);
  static int di[8];
  static int dj[8];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
};

#endif // __MOORE_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
