#include "olim.test.common.hpp"

#include <olim.hpp>

int main() {
  using olim = olim4_mp1;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();
}
