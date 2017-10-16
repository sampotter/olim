#ifndef __OLIM26_MP0_HPP__
#define __OLIM26_MP0_HPP__

#include "node_3d.hpp"
#include "olim26_rect.hpp"
#include "olim_rect_update_rules.hpp"
#include "speed_estimates.hpp"

using olim26_mp0 = olim26_rect<
	node_3d, olim_rect_update_rules, mp0_speed_estimate>;

#endif // __OLIM26_MP0_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
