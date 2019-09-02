#include "fac.h"

status_e fac_src_init(fac_src_s **f_ptr, int n, double const *x, double s)
{
  fac_src * f = (*f_ptr = new fac_src);

  f->n = n;

  f->x = new double[n];
  for (int i = 0; i < n; ++i) {
    f->x[i] = x[i];
  }

  f->s = s;

  return SUCCESS;
}

status_e fac_src_deinit(fac_src **f_ptr)
{
  fac_src * f = *f_ptr;

  delete f->x;

  delete *f_ptr;

  *f_ptr = nullptr;

  return SUCCESS;
}
