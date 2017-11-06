#ifndef __ABSTRACT_NODE_HPP__
#define __ABSTRACT_NODE_HPP__

#include <limits>

enum class state {valid, trial, far};

struct abstract_node {
  abstract_node() {}
  abstract_node(double value, state s): _value {value}, _state {s} {}
  inline double get_value() const { return _value; }
  inline void set_value(double value) { _value = value; }
  inline int get_heap_pos() const { return _heap_pos; }
  inline void set_heap_pos(int pos) { _heap_pos = pos; }
  inline bool is_valid() const { return _state == state::valid; }
  inline bool is_trial() const { return _state == state::trial; }
  inline bool is_far() const { return _state == state::far; }
  inline void set_valid() { _state = state::valid; }
  inline void set_trial() { _state = state::trial; }
  inline void set_far() { _state = state::far; }
protected:
  double _value {std::numeric_limits<double>::infinity()};
  state _state {state::far};
  int _heap_pos {-1};
};

#endif // __ABSTRACT_NODE_HPP__
