#include "olim4.hpp"
#include "olim.test.common.hpp"

int main() {
  trivial_case_works<olim4_rhr>();
  adjacent_update_works<olim4_rhr>();
}
