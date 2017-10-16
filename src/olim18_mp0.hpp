#ifndef __OLIM18_MP0_HPP__
#define __OLIM18_MP0_HPP__

#include "node_3d.hpp"
#include "olim18.hpp"
#include "olim_update_rules.hpp"
#include "speed_estimates.hpp"

using olim18_mp0 = olim18_rect<
  node_3d, olim_rect_update_rules, mp0_speed_estimate>;

#endif // __OLIM18_MP0_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
