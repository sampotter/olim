#ifndef __NODE_HPP__
#define __NODE_HPP__

#include <limits>

enum class state {valid, trial, far};

struct node {
  static node make_boundary_node(int i, int j, double value = 0.0);

  inline double get_value() const {
    return _value;
  }
  void set_value(double value);

  int get_i() const;
  void set_i(int i);
	
  int get_j() const;
  void set_j(int j);

  int get_heap_pos() const;
  void set_heap_pos(int pos);

  bool is_valid() const;
  bool is_trial() const;
  bool is_far() const;
  void set_valid();
  void set_trial();
  void set_far();
private:
  double _value {std::numeric_limits<double>::infinity()};
  state _state {state::far};
  int _i;
  int _j;
  int _heap_pos;
};

#endif // __NODE_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
