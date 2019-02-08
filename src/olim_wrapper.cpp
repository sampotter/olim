#include "olim_wrapper.h"

#include <string>
#include <variant>

#include "olim.hpp"
#include "olim3d.hpp"

struct olim_wrapper {
  olim4_rhr * olim;
};

status_e olim_wrapper_init(olim_wrapper_s ** w_ptr, olim_wrapper_params_s * p)
{
  olim_wrapper_s * w = (*w_ptr = new olim_wrapper_s);

  w->olim = new olim4_rhr {{p->dims}, p->h, no_slow_t {}};

  return SUCCESS;
}

status_e olim_wrapper_deinit(olim_wrapper_s ** w_ptr)
{
  olim_wrapper_s * w = *w_ptr;
  delete w->olim;
  
  delete *w_ptr;

  return SUCCESS;
}

status_e olim_wrapper_run(olim_wrapper_s * w)
{
  w->olim->run();

  return SUCCESS;
}

status_e olim_wrapper_add_bd(olim_wrapper_s * w, int * inds, double U)
{
  w->olim->add_boundary_node(vec2<int> {inds}, U);

  return SUCCESS;
}

status_e olim_wrapper_get_U_ptr(olim_wrapper_s * w, double ** U_ptr)
{
  *U_ptr = w->olim->get_U_ptr();

  return SUCCESS;
}

status_e olim_wrapper_get_s_ptr(olim_wrapper_s * w, double ** s_ptr)
{
  *s_ptr = w->olim->get_s_ptr();

  return SUCCESS;
}


status_e olim_wrapper_get_state_ptr(olim_wrapper_s * w, char ** state_ptr)
{
  *state_ptr = w->olim->get_state_ptr();

  return SUCCESS;
}
