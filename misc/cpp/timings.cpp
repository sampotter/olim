#include <chrono>
#include <iostream>
#include <vector>

#include <olim.hpp>
#include <olim3d.hpp>

using namespace std::chrono;

int main() {
  /**
   * Initialize problem sizes
   */
  std::vector<int> ns;
  for (int i = 3; i <= 8; ++i) {
    ns.push_back((1 << i) + 1);
  }
  
  for (auto n: ns) {
    double h = 2./(n - 1);
    int i0 = n/2;

    auto start = high_resolution_clock::now();

    olim26_rhr m {n, n, n, h, (speed_func_3d) s4, 1., 1., 1.};
    m.add_boundary_node(i0, i0, i0);
    m.run();

    auto stop = high_resolution_clock::now();

    auto t = duration<double>(stop - start).count();

    std::cout << n << ": " << t << "s" << std::endl;

  }
}
