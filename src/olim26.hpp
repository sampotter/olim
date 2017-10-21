#ifndef __OLIM26_HPP__
#define __OLIM26_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"
#include "olim_rect_update_rules.hpp"
#include "speed_estimates.hpp"

template <class node, class update_rules, class speed_estimate>
struct olim26_rect: public marcher_3d<node>, public update_rules,
                    public speed_estimate {
  using marcher_3d<node>::marcher_3d;
protected:
  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb);
  static int di[26];
  static int dj[26];
  static int dk[26];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
  virtual void update_impl(int i, int j, int k, double & T);
};

using olim26_rhr = olim26_rect<node_3d, olim_rect_update_rules, rhr_speed_estimate>;
using olim26_mp0 = olim26_rect<node_3d, olim_rect_update_rules, mp0_speed_estimate>;

#include "olim26.impl.hpp"

#endif // __OLIM26_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
