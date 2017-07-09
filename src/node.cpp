#include "node.hpp"

node node::make_boundary_node(int i, int j, double value) {
  node n;
  n._value = value;
  n._state = state::valid;
  n._i = i;
  n._j = j;
  return n;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
