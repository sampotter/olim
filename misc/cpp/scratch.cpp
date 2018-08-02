#include <olim.hpp>
#include <olim3d.hpp>
#include <src/config.hpp>
#include <string>

int main(int argc, char * argv[]) {
  if (argc < 3) {
    std::exit(0);
  }

  std::string marcher_name {argv[1]};
  int n = atoi(argv[2]);

  if (marcher_name == "olim3d_hu_mp1") {
    olim3d_hu_mp1 m {n, n, n, 2./(n-1), (speed_func_3d) s1, 1., 1., 1.};
    m.add_boundary_node(n/2, n/2, n/2);
    m.run();
#if COLLECT_STATS
    m.dump_stats();
#endif
  } else if (marcher_name == "olim26_mp1") {
    olim26_mp1 m {n, n, n, 2./(n-1), (speed_func_3d) s1, 1., 1., 1.};
    m.add_boundary_node(n/2, n/2, n/2);
    m.run();
#if COLLECT_STATS
    m.dump_stats();
#endif
  }
}
