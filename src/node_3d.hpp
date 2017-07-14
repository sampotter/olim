#ifndef __NODE_3D_HPP__
#define __NODE_3D_HPP__

#include <limits>

enum class state {valid, trial, far};

struct node_3d {
  static node_3d make_boundary_node(int i, int j, int k, double value = 0.0);

  inline double get_value() const {
    return _value;
  }

  inline void set_value(double value) {
    _value = value;
  }

  inline int get_i() const {
    return _i;
  }

  inline void set_i(int i) {
    _i = i;
  }
	
  inline int get_j() const {
    return _j;
  }

  inline void set_j(int j) {
    _j = j;
  }

  inline int get_k() const {
    return _k;
  }

  inline void set_k(int k) {
    _k = k;
  }

  inline int get_heap_pos() const {
    return _heap_pos;
  }

  inline void set_heap_pos(int pos) {
    _heap_pos = pos;
  }

  inline bool is_valid() const {
    return _state == state::valid;
  }

  inline bool is_trial() const {
    return _state == state::trial;
  }

  inline bool is_far() const {
    return _state == state::far;
  }

  inline void set_valid() {
    _state = state::valid;
  }

  inline void set_trial() {
    _state = state::trial;
  }

  inline void set_far() {
    _state = state::far;
  }
private:
  double _value {std::numeric_limits<double>::infinity()};
  state _state {state::far};
  int _i;
  int _j;
  int _k;
  int _heap_pos;
};

#endif // __NODE_3D_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
