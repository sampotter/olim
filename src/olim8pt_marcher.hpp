#ifndef __OLIM8PT_MARCHER_HPP__
#define __OLIM8PT_MARCHER_HPP__

#include "moore_marcher.hpp"

struct olim8pt_marcher: public moore_marcher
{
  using moore_marcher::moore_marcher;
private:
  virtual void update_node_value_impl(size_t i, size_t j);
  double get_alpha(double u0, double u1, double s) const;
  double solve2pt_adjacent(double u0, double u1, double sh) const;
  double solve2pt_diagonal(double u0, double u1, double sh) const;
};

#endif // __OLIM8PT_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
