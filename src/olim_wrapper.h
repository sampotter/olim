#ifndef OLIM_WRAPPER_h
#define OLIM_WRAPPER_h

#ifdef __cplusplus
extern "C" {
#endif

#include "type.h"

typedef enum status {
  SUCCESS
} status_e;
  
typedef struct fac_src_wrapper_params {
  int ndims;
  double *coords;
  double s;
} fac_src_wrapper_params_s;

typedef struct fac_src_wrapper fac_src_wrapper_s;

status_e fac_src_wrapper_init(fac_src_wrapper_s **w_ptr, fac_src_wrapper_params_s *p);
status_e fac_src_wrapper_deinit(fac_src_wrapper_s **w_ptr);

typedef struct olim_wrapper_params {
  enum neighborhood nb;
  enum cost_func F;
  double h;
  int ndims;
  int *dims;
} olim_wrapper_params_s;
  
typedef struct olim_wrapper olim_wrapper_s;

status_e olim_wrapper_init(olim_wrapper_s **w_ptr, olim_wrapper_params_s *p);
status_e olim_wrapper_deinit(olim_wrapper_s **w_ptr);
status_e olim_wrapper_solve(olim_wrapper_s *w);
status_e olim_wrapper_step(olim_wrapper_s *w, int *lin);
status_e olim_wrapper_peek(olim_wrapper_s *w, double *value, int *lin, bool *empty);
status_e olim_wrapper_adjust(olim_wrapper_s *w, int *inds, double U);
status_e olim_wrapper_add_src(olim_wrapper_s *w, int *inds, double U);
status_e olim_wrapper_add_bd(olim_wrapper_s *w, int *inds);
status_e olim_wrapper_add_free(olim_wrapper_s *w, int *inds);
status_e olim_wrapper_set_fac_src(olim_wrapper_s *w, int *inds, fac_src_wrapper_s *fs);
status_e olim_wrapper_get_U_ptr(olim_wrapper_s *w, double **U_ptr);
status_e olim_wrapper_get_s_ptr(olim_wrapper_s *w, double **s_ptr);
status_e olim_wrapper_get_state_ptr(olim_wrapper_s *w, char **state_ptr);

#ifdef __cplusplus
}
#endif

#endif // OLIM_WRAPPER_h
