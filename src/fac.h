#ifndef __FAC_H__
#define __FAC_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "type.h"

  struct fac_src
  {
    int n;

    /**
     * Location of factored source in "index space"---may not be
     * grid-aligned, hence doubles.
     */
    double * inds;

    double s;
  };

  typedef struct fac_src fac_src_s;

  /**
   * Initialize a factored source:
   *
   * - n: the dimension of the problem (n = 2, 3 for now)
   * - inds: points to n doubles corresponding to the location of the
   *   source in "index space"
   * - s: the value of the slowness at the factored source
   * - padding: a possible padding value to offset `inds` by if the
   *   underlying marcher uses padding for speed (e.g., in the case of
   *   the eikonal solver)
   *
   * TODO: later, we may want to replace `s` with a void* to user
   * defined data if a double isn't enough to specify the boundary
   * conditions at the factored source.
   */
  status_e fac_src_init(fac_src_s **f_ptr, int n, double const *inds, double s,
                        int const *padding);

  status_e fac_src_deinit(fac_src **f_ptr);

#ifdef __cplusplus
}
#endif

#endif // __FAC_H__
