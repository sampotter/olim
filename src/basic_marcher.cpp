#include "basic_marcher.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

void basic_marcher::update_node_value_impl(int i, int j, double & T) {
  abstract_node * nb[4] = {nullptr, nullptr, nullptr, nullptr};
  get_valid_neighbors(i, j, nb);
  double sh = get_h()*S(i, j);

  double T1 = std::min(
    nb[0] ? nb[0]->get_value() : std::numeric_limits<double>::infinity(),
    nb[2] ? nb[2]->get_value() : std::numeric_limits<double>::infinity());

  double T2 = std::min(
    nb[1] ? nb[1]->get_value() : std::numeric_limits<double>::infinity(),
    nb[3] ? nb[3]->get_value() : std::numeric_limits<double>::infinity());

  bool T1_inf = std::isinf(T1), T2_inf = std::isinf(T2);

  if (!T1_inf && !T2_inf) {
    double diff = T1 - T2, disc = 2*sh*sh - diff*diff;
    T = disc > 0 ? std::min(T, (T1 + T2 + std::sqrt(disc))/2) : T;
  } else if (std::isinf(T1)) {
    T = std::min(T, T2 + sh);
  } else if (std::isinf(T2)) {
    T = std::min(T, T1 + sh);
  } else {
    assert(false);
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
