#ifndef OLIM_WRAPPER_h
#define OLIM_WRAPPER_h

#ifdef __cplusplus
extern "C" {
#endif

#include "type.h"

typedef enum status {
  SUCCESS
} status_e;
  
typedef struct olim_wrapper_params {
  enum neighborhood nb;
  enum cost_func F;
  double h;
  int ndims;
  int * dims;
} olim_wrapper_params_s;
  
typedef struct olim_wrapper olim_wrapper_s;

status_e olim_wrapper_init(olim_wrapper_s ** w_ptr, olim_wrapper_params_s * p);
status_e olim_wrapper_deinit(olim_wrapper_s ** w_ptr);
status_e olim_wrapper_run(olim_wrapper_s * w);
status_e olim_wrapper_add_bd(olim_wrapper_s * w, int * inds, double U);
status_e olim_wrapper_get_U_ptr(olim_wrapper_s * w, double ** U_ptr);
status_e olim_wrapper_get_s_ptr(olim_wrapper_s * w, double ** s_ptr);
  status_e olim_wrapper_get_state_ptr(olim_wrapper_s * w, char ** state_ptr);

#ifdef __cplusplus
}
#endif

#endif // OLIM_WRAPPER_h
