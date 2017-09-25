#ifndef __OLIM26_RHR_HPP__
#define __OLIM26_RHR_HPP__

#include "conics.hpp"
#include "node_3d.hpp"
#include "olim26.hpp"
#include "olim_update_rules.hpp"

using olim26_rhr_arma = olim26<
  node_3d, olim3d_rhr_update_rules<arma_rootfinder>>;

#endif // __OLIM26_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
