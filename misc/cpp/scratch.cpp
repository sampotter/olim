#include <olim.hpp>
#include <olim3d.hpp>
#include <src/config.hpp>
#include <string>

constexpr int n = 65;

int main() {

  int i0 = n/2;
  double h = 2./(n-1);

  double * S = new double[n*n*n];
  for (int i = 0; i < n*n*n; ++i) S[i] = 1;
  
  olim26_mp0 o {n, n, n, h, S};
  // olim3d_hu_mp0 o {n, n, n, h, S};

  // for (int i = 0; i < n; ++i) {
  //   for (int j = 0; j < n; ++j) {
  //     for (int k = 0; k < n; ++k) {
  //       o.set_node_fac_parent(i, j, k, i0, i0, i0);
  //     }
  //   }
  // }

  o.add_boundary_node(i0, i0, i0);
  o.run();
}
