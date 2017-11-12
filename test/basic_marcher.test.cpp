#include "basic_marcher.hpp"
#include "olim.test.common.hpp"

int main() {
  using olim = basic_marcher;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();
  error_is_monotonic<olim>();
}
