#include "olim4.hpp"
#include "olim.test.common.hpp"

int main() {
  using olim = olim4_rhr;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();
  error_is_monotonic<olim>();
}
