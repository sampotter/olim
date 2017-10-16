#ifndef __OLIM18_HPP__
#define __OLIM18_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"
#include "olim_update_rules.hpp"
#include "speed_estimates.hpp"

template <class node, class update_rules, class speed_estimate>
struct olim18_rect: public marcher_3d<node>, public update_rules,
                    public speed_estimate {
  using marcher_3d<node>::marcher_3d;
protected:
  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb);
  static int di[18];
  static int dj[18];
  static int dk[18];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
  virtual void update_impl(int i, int j, int k, double & T);
  void do_tri12_updates(
    abstract_node const * const * nb, int const * dirs,
    double const * s_, double s, double h, double & T) const;
  void do_tri22_updates(
    abstract_node const * const * nb, int const * dirs,
    double const * s_, double s, double h, double & T) const;
};

using olim18_mp0 = olim18_rect<node_3d, olim_rect_update_rules, mp0_speed_estimate>;
using olim18_rhr = olim18_rect<node_3d, olim_rect_update_rules, rhr_speed_estimate>;

#include "olim18.impl.hpp"

#endif // __OLIM18_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
