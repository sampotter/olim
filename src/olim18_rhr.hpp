#ifndef __OLIM18_RHR_HPP__
#define __OLIM18_RHR_HPP__

#include "conics.hpp"
#include "node_3d.hpp"
#include "olim18.hpp"
#include "olim_update_rules.hpp"

using olim18_rhr_arma = olim18<
  node_3d, olim3d_rhr_update_rules<arma_rootfinder>>;

#endif // __OLIM18_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
