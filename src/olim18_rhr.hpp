#ifndef __OLIM18_RHR_HPP__
#define __OLIM18_RHR_HPP__

#include "node_3d.hpp"
#include "olim18_rect.hpp"
#include "olim_update_rules.hpp"
#include "speed_estimates.hpp"

using olim18_rhr = olim18_rect<
  node_3d, olim_rect_update_rules, rhr_speed_estimate>;

#endif // __OLIM18_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
