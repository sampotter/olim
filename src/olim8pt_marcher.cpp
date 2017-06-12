#include "olim8pt_marcher.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

void olim8pt_marcher::update_node_value_impl(size_t i, size_t j) {
  node* nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  node *x = 0x0, *x0 = 0x0, *x1 = 0x0;
  get_valid_neighbors(i, j, nb);
  double T = std::numeric_limits<double>::infinity();
  double sh = get_h()/F(i, j);

  // TODO: cache alphas?

  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 6) % 8] && !nb[(k + 7) % 8] &&
        !nb[(k + 1) % 8] && !nb[(k + 2) % 8]) {
      T = std::min(T, x0->get_value() + sh);
    }
  }
  for (int k = 1; k < 8; k += 2) {
    if ((x0 = nb[k]) && !nb[(k + 7) % 8] && !nb[(k + 1) % 8]) {
      T = std::min(T, x0->get_value() + sh*sqrt(2));
    }
  }

  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && (x1 = nb[k + 1])) {
      T = std::min(T, solve2pt_diagonal(x0->get_value(), x1->get_value(), sh));
    }
  }
  for (int k = 1; k < 8; k += 2) {
    if ((x0 = nb[(k + 1) % 8]) && (x1 = nb[k])) {
      T = std::min(T, solve2pt_diagonal(x0->get_value(), x1->get_value(), sh));
    }
  }

  for (int k = 0; k < 8; k += 2) {
    if ((x0 = nb[k]) && (x1 = nb[(k + 2) % 8])) {
      T = std::min(T, solve2pt_adjacent(x0->get_value(), x1->get_value(), sh));
    }
  }

  // TODO: move this out one level

  x = &this->operator()(i, j);
  assert(x->is_trial());
  if (T <= x->get_value()) {
    x->set_value(T);
    adjust_heap_entry(x);
  }
}

double olim8pt_marcher::get_alpha(double u0, double u1, double sh) const {
  return std::fabs(u0 - u1)/sh;
}

double olim8pt_marcher::solve2pt_adjacent(double u0, double u1, double sh)
  const
{
  assert(u0 >= 0);
  assert(u1 >= 0);
  assert(!std::isinf(u0));
  assert(!std::isinf(u1));
  double alpha = get_alpha(u0, u1, sh);
  if (alpha > std::sqrt(2)/2) {
    return std::numeric_limits<double>::infinity();
  }
  double rad = sqrt(3)*alpha/sqrt(alpha*alpha + 2);
  double lam1 = (1 + rad)/2;
  double lam2 = (1 - rad)/2;
  assert(std::fabs(1 - lam1 - lam2) < 1e-15 ||
         (0 < lam1 && lam1 < 1) != (0 < lam2 && lam2 < 1));
  double lam = 0 < lam1 && lam1 < 1 ? lam1 : lam2;
  return (1 - lam)*u0 + lam*u1 + sh*sqrt(2*lam*(1 - lam) + 1);
}

double olim8pt_marcher::solve2pt_diagonal(double u0, double u1, double sh)
  const
{
  assert(u0 >= 0);
  assert(u1 >= 0);
  assert(!std::isinf(u0));
  assert(!std::isinf(u1));
  double alpha = get_alpha(u0, u1, sh);
  if (alpha > std::sqrt(2)/2) {
    return std::numeric_limits<double>::infinity();
  }
  double lam = alpha/(1 - alpha*alpha);
  assert(0 <= lam && lam <= 1);
  return (1 - lam)*u0 + lam*u1 + sh*sqrt(lam*lam + 1);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
