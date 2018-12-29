#ifndef __ABSTRACT_NODE_HPP__
#define __ABSTRACT_NODE_HPP__

#include <src/config.hpp>

#include <array>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>

#include "common.hpp"

enum class state {valid, trial, far};

inline
std::string
to_string(state const & s)
{
  if (s == state::valid) return "valid";
  else if (s == state::trial) return "trial";
  else return "far";
}

struct abstract_node {
  abstract_node(): _fac_parent {nullptr} {}
  abstract_node(double value, state s):
    _value {value}, _state {s}, _fac_parent {nullptr} {}
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
  inline bool has_fac_parent() const { return _fac_parent != nullptr; }
  inline abstract_node * get_fac_parent() const { return _fac_parent; }
  inline void set_fac_parent(abstract_node * parent) { _fac_parent = parent; }
#if TRACK_PARENTS
  inline bool has_parents() const {
    return _parents[0] != nullptr && _parents[1] != nullptr &&
      _parents[2] != nullptr;
  }
  inline std::array<abstract_node *, 3> get_parents() const { return _parents; }
  inline void set_parents(std::array<abstract_node *, 3> const & parents) {
    _parents = parents;
  }
#endif
#if NODE_MONITORING
  inline void start_monitoring() { _monitor_node = true; }
  inline void stop_monitoring() { _monitor_node = false; }
  inline bool monitoring_node() const { return _monitor_node; }
#endif
EIKONAL_PROTECTED:
  double _value {std::numeric_limits<double>::infinity()};
  state _state {state::far};
  int _heap_pos {-1};
  abstract_node * _fac_parent;
#if TRACK_PARENTS
  std::array<abstract_node *, 3> _parents {{nullptr, nullptr, nullptr}};
#endif
#if NODE_MONITORING
  bool _monitor_node {false};
#endif
};

inline
std::string
to_string(abstract_node const & n)
{
  std::ostringstream os;
  os.precision(std::numeric_limits<double>::max_digits10);
  os << "value = " << n.get_value() << ", "
     << "state = " << to_string(n.get_state());
  return os.str();
}

inline
std::ostream &
operator<<(std::ostream & o, abstract_node const & n)
{
  return o << "abstract_node {" << to_string(n) << "}";
}

#endif // __ABSTRACT_NODE_HPP__
