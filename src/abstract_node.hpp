#ifndef __ABSTRACT_NODE_HPP__
#define __ABSTRACT_NODE_HPP__

#include <limits>

#include "common.macros.hpp"

enum class state {valid, trial, far};

struct abstract_node {
  abstract_node(): _parent {nullptr} {}
  abstract_node(double value, state s):
    _value {value}, _state {s}, _parent {nullptr} {}
  inline double get_value() const { return _value; }
  inline void set_value(double value) { _value = value; }
  inline int get_heap_pos() const { return _heap_pos; }
  inline void set_heap_pos(int pos) { _heap_pos = pos; }
  inline state get_state() const { return _state; }
  inline void set_state(state s) { _state = s; }
  inline bool is_valid() const { return _state == state::valid; }
  inline bool is_trial() const { return _state == state::trial; }
  inline bool is_far() const { return _state == state::far; }
  inline void set_valid() { _state = state::valid; }
  inline void set_trial() { _state = state::trial; }
  inline void set_far() { _state = state::far; }
  inline bool has_parent() const { return _parent != nullptr; }
  inline abstract_node * get_parent() const { return _parent; }
  inline void set_parent(abstract_node * parent) { _parent = parent; }
  inline double get_sh_fac() const { return _sh_fac; }
  inline void set_sh_fac(double sh_fac) { _sh_fac = sh_fac; }
  inline virtual double get_T() const = 0;
EIKONAL_PROTECTED:
  double _value {std::numeric_limits<double>::infinity()};
  state _state {state::far};
  int _heap_pos {-1};
  union {
    abstract_node * _parent;
    double _sh_fac;
  };
};

#endif // __ABSTRACT_NODE_HPP__
