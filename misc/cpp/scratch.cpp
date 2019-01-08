#include <olim.hpp>
#include <olim3d.hpp>
#include <src/config.hpp>
#include <string>

constexpr int n = 33;

int main() {

  int i0 = n/2;
  double h = 2./(n-1);

  double * S = new double[n*n*n];
  // for (int i = 0; i < n*n*n; ++i) S[i] = 1;
  double x, y, z;
  for (int i = 0; i < n; ++i) {
    x = h*static_cast<double>(i) - 1;
    for (int j = 0; j < n; ++j) {
      y = h*static_cast<double>(j) - 1;
      for (int k = 0; k < n; ++k) {
        z = h*static_cast<double>(k) - 1;
        S[n*(n*k + j) + i] = s2(x, y, z);
      }
    }
  }
  
  // olim18_mp1 o {n, n, n, h, S};
  // olim26_mp0 o {n, n, n, h, S};
  olim3d_hu_mp0 o {n, n, n, h, S};

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
