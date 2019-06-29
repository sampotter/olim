#ifndef QUASIPOT_H
#define QUASIPOT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct quasipot_params {
  int ndims;
  int *dims;
  double h;
  void *b;
  int K;
};

typedef struct quasipot quasipot_s;

status_e quasipot_init(quasipot_s **handle, quasipot_params_s *p);
status_e quasipot_deinit(quasipot_s **handle);
status_e quasipot_solve(quasipot_s *q);
status_e quasipot_step(quasipot_s *q, int *lin);
status_e quasipot_peek(quasipot_s *q, double *value, int *lin, bool *empty);
status_e quasipot_adjust(quasipot_s *q, int *inds, double U);
status_e quasipot_add_src(quasipot_s *q, int *inds, double U);
status_e quasipot_add_bd(quasipot_s *q, int *inds);
status_e quasipot_add_free(quasipot_s *q, int *inds);
status_e quasipot_get_U_ptr(quasipot_s *q, double **U_ptr);
status_e quasipot_get_b_ptr(quasipot_s *q, double **b_ptr);
status_e quasipot_get_state_ptr(quasipot_s *q, char **state_ptr);

#ifdef __cplusplus
}
#endif

#endif // QUASIPOT_WRAPPER_H
