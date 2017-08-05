#include "basic_marcher_3d.hpp"

int main() {
  basic_marcher_3d m {3, 3, 3};
  m.add_boundary_node(1, 1, 1);
  m.run();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
