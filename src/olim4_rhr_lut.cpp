#include "olim4_rhr_lut.hpp"

#include <cassert>

#include "olim_util.hpp"

constexpr int NB_W = 1;
constexpr int NB_S = 2;
constexpr int NB_E = 4;
constexpr int NB_N = 8;

void olim4_rhr_lut::update_node_value_impl(int i, int j, double & T) {
  node * nb[4] = {0x0, 0x0, 0x0, 0x0}; // NESW
  get_valid_neighbors(i, j, (abstract_node **) nb);
  double s = S(i, j), h = get_h(), sh = s*h, tmp;

  switch ((nb[0] ? NB_N : 0) | (nb[1] ? NB_E : 0) | (nb[2] ? NB_S : 0) |
          (nb[3] ? NB_W : 0)) {
  /* one neighbor */
  case NB_N:
    T = nb[0]->get_value() + sh;
    break;
  case NB_E:
    T = nb[1]->get_value() + sh;
    break;
  case NB_S:
    T = nb[2]->get_value() + sh;
    break;
  case NB_W:
    T = nb[3]->get_value() + sh;
    break;
  /* two neighbors */
  case NB_N | NB_E:
    T = rhr_adj(nb[0]->get_value(), nb[1]->get_value(), s, h);
    break;
  case NB_N | NB_S:
    T = std::min(nb[0]->get_value() + sh, nb[2]->get_value() + sh);
    break;
  case NB_N | NB_W:
    T = rhr_adj(nb[0]->get_value(), nb[3]->get_value(), s, h);
    break;
  case NB_E | NB_S:
    T = rhr_adj(nb[1]->get_value(), nb[2]->get_value(), s, h);
    break;
  case NB_E | NB_W:
    T = std::min(nb[1]->get_value() + sh, nb[3]->get_value() + sh);
    break;
  case NB_S | NB_W:
    T = rhr_adj(nb[2]->get_value(), nb[3]->get_value(), s, h);
    break;
  /* three neighbors */
  case NB_N | NB_E | NB_S:
    tmp = nb[1]->get_value();
    T = std::min(
      rhr_adj(nb[0]->get_value(), tmp, s, h),
      rhr_adj(tmp, nb[2]->get_value(), s, h));
    break;
  case NB_N | NB_E | NB_W:
    tmp = nb[0]->get_value();
    T = std::min(
      rhr_adj(tmp, nb[1]->get_value(), s, h),
      rhr_adj(nb[3]->get_value(), tmp, s, h));
    break;
  case NB_N | NB_S | NB_W:
    tmp = nb[3]->get_value();
    T = std::min(
      rhr_adj(nb[2]->get_value(), tmp, s, h),
      rhr_adj(tmp, nb[0]->get_value(), s, h));
    break;
  case NB_E | NB_S | NB_W:
    tmp = nb[2]->get_value();
    T = std::min(
      rhr_adj(nb[1]->get_value(), tmp, s, h),
      rhr_adj(tmp, nb[3]->get_value(), s, h));
    break;
  /* four neighbors */
  case NB_N | NB_E | NB_S | NB_W:
    tmp = nb[1]->get_value();
    T = std::min(
      std::min(
        rhr_adj(nb[0]->get_value(), tmp, s, h),
        rhr_adj(tmp, nb[2]->get_value(), s, h)),
      std::min(
        rhr_adj(nb[2]->get_value(), nb[3]->get_value(), s, h),
        rhr_adj(nb[3]->get_value(), nb[0]->get_value(), s, h)));
    break;
  default:
    assert(false);
    break;
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
