#ifndef __FAC_H__
#define __FAC_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "type.h"

  struct fac_src
  {
    int n;
    double * x;
    double s;
  };

  typedef struct fac_src fac_src_s;

  status_e fac_src_init(fac_src_s **f_ptr, int n, double const *x, double s);
  status_e fac_src_deinit(fac_src **f_ptr);

#ifdef __cplusplus
}
#endif

#endif // __FAC_H__
