#include <basic_marcher.hpp>
#include <basic_marcher_3d.hpp>
#include <olim.hpp>
#include <olim3d.hpp>

#include <iostream>
#include <string>

template <class marcher_3d>
void run_marcher_3d(int n) {
  double h = 2./(n - 1);
  int i = n/2;
  marcher_3d m {n, n, n, h, (speed_func_3d) default_speed_func, 1., 1., 1.};
  m.add_boundary_node(i, i, i);
  m.run();
  m.dump_stats();
}

int main(int argc, char * argv[]) {
  if (argc != 3) {
    std::cout << "usage: " << argv[0] << " marcher n" << std::endl;
    std::exit(1);
  }

  std::string marcher_name {argv[1]};
  int n = std::stoi(argv[2]);

  if (marcher_name == "olim6_mp0") run_marcher_3d<olim6_mp0>(n);
  if (marcher_name == "olim6_mp1") run_marcher_3d<olim6_mp1>(n);
  if (marcher_name == "olim6_rhr") run_marcher_3d<olim6_rhr>(n);

  if (marcher_name == "olim18_mp0") run_marcher_3d<olim18_mp0>(n);
  if (marcher_name == "olim18_mp1") run_marcher_3d<olim18_mp1>(n);
  if (marcher_name == "olim18_rhr") run_marcher_3d<olim18_rhr>(n);

  if (marcher_name == "olim26_mp0") run_marcher_3d<olim26_mp0>(n);
  if (marcher_name == "olim26_mp1") run_marcher_3d<olim26_mp1>(n);
  if (marcher_name == "olim26_rhr") run_marcher_3d<olim26_rhr>(n);

  if (marcher_name == "olim3d_hu_mp0") run_marcher_3d<olim3d_hu_mp0>(n);
  if (marcher_name == "olim3d_hu_mp1") run_marcher_3d<olim3d_hu_mp1>(n);
  if (marcher_name == "olim3d_hu_rhr") run_marcher_3d<olim3d_hu_rhr>(n);
}
