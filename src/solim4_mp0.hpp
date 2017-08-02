#ifndef __SOLIM4_MP0_HPP__
#define __SOLIM4_MP0_HPP__

#include "smart_neumann_marcher.hpp"

struct solim4_mp0: public smart_neumann_marcher {
  using smart_neumann_marcher::smart_neumann_marcher;
private:
  void update_impl(int i, int j, double & T) override final;
};

#endif // __SOLIM4_MP0_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
