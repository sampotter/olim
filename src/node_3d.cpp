#include "node_3d.hpp"

node_3d node_3d::make_boundary_node(int i, int j, int k, double value) {
  node_3d n;
  n._value = value;
  n._state = state::valid;
  n._i = i;
  n._j = j;
  n._k = k;
  return n;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
