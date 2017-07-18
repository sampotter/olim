#ifndef __SOLIM4_MP0_HPP__
#define __SOLIM4_MP0_HPP__

#include "neumann_marcher.hpp"
#include "smart_node.hpp"

struct solim4_mp0: public neumann_marcher<smart_node> {
  using neumann_marcher::neumann_marcher;
private:
  void update_node_value_impl(int i, int j, double & T) override final;
};

#endif // __SOLIM4_MP0_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
