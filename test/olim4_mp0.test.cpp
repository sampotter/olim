#include "olim4.hpp"
#include "olim.test.common.hpp"

int main() {
  using olim = olim4_mp0;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();
}