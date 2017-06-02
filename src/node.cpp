#include "node.hpp"

node node::make_boundary_node(size_t i, size_t j) {
  node n;
  n._value = 0;
  n._state = state::valid;
  n._i = i;
  n._j = j;
  return n;
}

double node::get_value() const {
  return _value;
}

void node::set_value(double value) {
  _value = value;
}

size_t node::get_i() const {
  return _i;
}

void node::set_i(size_t i) {
  _i = i;
}

size_t node::get_j() const {
  return _j;
}

void node::set_j(size_t j) {
  _j = j;
}

size_t node::get_heap_pos() const {
  return _heap_pos;
}

void node::set_heap_pos(size_t pos) {
  _heap_pos = pos;
}

bool node::is_valid() const {
  return _state == state::valid;
}

bool node::is_trial() const {
  return _state == state::trial;
}

bool node::is_far() const {
  return _state == state::far;
}

void node::set_valid() {
  _state = state::valid;
}

void node::set_trial() {
  _state = state::trial;
}

void node::set_far() {
  _state = state::far;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
