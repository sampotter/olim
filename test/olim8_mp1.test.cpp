#include "olim8.hpp"
#include "olim.test.common.hpp"

int main() {
  using olim = olim8_mp1;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();

  int nmin = 5, nmax = 31;
  for (size_t i = 0; i < speed_funcs.size(); ++i) {
    auto s = speed_funcs[i];
    auto f = speed_func_solns[i];
    error_is_monotonic<olim>(nmin, nmax, s, f);
  }
}
