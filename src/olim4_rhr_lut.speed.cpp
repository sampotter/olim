#include "olim4_rhr_lut.hpp"
#include "speed_funcs.hpp"

int main() {
  olim4_rhr_lut m {1025, 1025, 1./512., default_speed_func, 1., 1.};
  m.add_boundary_node(512, 512);
  m.run();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
