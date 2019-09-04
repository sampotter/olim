#include "fac.h"

status_e fac_src_init(fac_src_s **f_ptr, int n, double const *inds, double s,
                      int const *padding)
{
  fac_src * f = (*f_ptr = new fac_src);

  f->n = n;

  f->inds = new double[n];
  for (int i = 0; i < n; ++i) {
    f->inds[i] = inds[i];
  }

  if (padding != nullptr) {
    for (int i = 0; i < n; ++i) {
      f->inds[i] += *padding;
    }
  }

  f->s = s;

  return SUCCESS;
}

status_e fac_src_deinit(fac_src **f_ptr)
{
  fac_src * f = *f_ptr;

  delete f->inds;

  delete *f_ptr;

  *f_ptr = nullptr;

  return SUCCESS;
}
