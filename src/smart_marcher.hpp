#ifndef __SMART_MARCHER_HPP__
#define __SMART_MARCHER_HPP__

#include "marcher.hpp"
#include "smart_node.hpp"

struct smart_marcher: public marcher {
  smart_marcher(int height, int width, double h = 1,
                speed_func S = default_speed_func,
                double x0 = 0.0, double y0 = 0.0);
  smart_marcher(int height, int width, double h, double * S_values);
  
  void add_boundary_node(int i, int j, double value = 0.0) override final;
  double get_value(int i, int j) const override final;
  
protected:
  smart_node & operator()(int i, int j);
  smart_node const & operator()(int i, int j) const;
  void stage_neighbor(int i, int j);
  void update_node_value(int i, int j);
  bool is_valid(int i, int j) const;

private:
  void init();
  
  smart_node * _nodes;
};

#endif // __SMART_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
