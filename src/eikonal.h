#ifndef EIKONAL_H
#define EIKONAL_H

#ifdef __cplusplus
extern "C" {
#endif

#include "fac.h"
#include "type.h"

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
status_e eikonal_solve(eikonal_wkspc_s *w);
status_e eikonal_step(eikonal_wkspc_s *w, int *lin);
status_e eikonal_peek(eikonal_wkspc_s *w, double *value, int *lin, bool *empty);
status_e eikonal_adjust(eikonal_wkspc_s *w, int *inds, double U);
status_e eikonal_add_src(eikonal_wkspc_s *w, int *inds, double U);
status_e eikonal_add_bd(eikonal_wkspc_s *w, int *inds);
status_e eikonal_add_free(eikonal_wkspc_s *w, int *inds);
status_e eikonal_factor(eikonal_wkspc_s *w, int *inds, fac_src_s *fs);
status_e eikonal_get_U_ptr(eikonal_wkspc_s *w, double **U_ptr);
status_e eikonal_get_s_ptr(eikonal_wkspc_s *w, double **s_ptr);
status_e eikonal_get_state_ptr(eikonal_wkspc_s *w, char **state_ptr);

#ifdef __cplusplus
}
#endif

#endif // EIKONAL_H
