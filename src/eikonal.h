#ifndef EIKONAL_H
#define EIKONAL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "type.h"

typedef struct fac_src_wrapper_params {
  int ndims;
  double *coords;
  double s;
} fac_src_wrapper_params_s;

typedef struct fac_src_wrapper fac_src_wrapper_s;

status_e fac_src_wrapper_init(fac_src_wrapper_s **w_ptr, fac_src_wrapper_params_s *p);
status_e fac_src_wrapper_deinit(fac_src_wrapper_s **w_ptr);

typedef struct eikonal_params {
  enum neighborhood nb;
  enum cost_func F;
  double h;
  int ndims;
  int *dims;
} eikonal_params_s;

typedef struct eikonal_wkspc eikonal_wkspc_s;

status_e eikonal_init(eikonal_wkspc_s **w_ptr, eikonal_params_s *p);
status_e eikonal_deinit(eikonal_wkspc_s **w_ptr);
status_e eikonal_wkspc_solve(eikonal_wkspc_s *w);
status_e eikonal_wkspc_step(eikonal_wkspc_s *w, int *lin);
status_e eikonal_peek(eikonal_wkspc_s *w, double *value, int *lin, bool *empty);
status_e eikonal_adjust(eikonal_wkspc_s *w, int *inds, double U);
status_e eikonal_add_src(eikonal_wkspc_s *w, int *inds, double U);
status_e eikonal_add_bd(eikonal_wkspc_s *w, int *inds);
status_e eikonal_add_free(eikonal_wkspc_s *w, int *inds);
status_e eikonal_wkspc_set_fac_src(eikonal_wkspc_s *w, int *inds, fac_src_wrapper_s *fs);
status_e eikonal_get_U_ptr(eikonal_wkspc_s *w, double **U_ptr);
status_e eikonal_get_s_ptr(eikonal_wkspc_s *w, double **s_ptr);
status_e eikonal_get_state_ptr(eikonal_wkspc_s *w, char **state_ptr);

#ifdef __cplusplus
}
#endif

#endif // EIKONAL_H
