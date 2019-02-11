#include "olim_wrapper.h"

#include "olim.hpp"
#include "olim3d.hpp"

// TODO: this is pretty awful, but at least it works. Definitely need
// to improve on this somehow. I tried std::variant, but ran into
// problems with it. Not totally sure what the cause of them was.

struct null_olim {
  void run() {}
  void add_src(int *, double) {}
  void add_bd(int *, double) {}
  double * get_U_ptr() const { return nullptr; }
  double * get_s_ptr() const { return nullptr; }
  char * get_state_ptr() const { return nullptr; }
};

union olim_variant {
  ~olim_variant() {}
  null_olim _null_olim;
  olim4_mp0 _olim4_mp0;
  olim4_mp1 _olim4_mp1;
  olim4_rhr _olim4_rhr;
  olim8_mp0 _olim8_mp0;
  olim8_mp1 _olim8_mp1;
  olim8_rhr _olim8_rhr;
  olim6_mp0 _olim6_mp0;
  olim6_mp1 _olim6_mp1;
  olim6_rhr _olim6_rhr;
  olim18_mp0 _olim18_mp0;
  olim18_mp1 _olim18_mp1;
  olim18_rhr _olim18_rhr;
  olim26_mp0 _olim26_mp0;
  olim26_mp1 _olim26_mp1;
  olim26_rhr _olim26_rhr;
  olim3d_mp0 _olim3d_mp0;
  olim3d_mp1 _olim3d_mp1;
  olim3d_rhr _olim3d_rhr;
};

struct olim_wrapper
{
  neighborhood nb;
  cost_func F;
  olim_variant _olim;

  template <class olim_t>
  olim_t & olim();

  template <>
  olim4_mp0 & olim<olim4_mp0>() {
    return _olim._olim4_mp0;
  }

  template <>
  olim4_mp1 & olim<olim4_mp1>() {
    return _olim._olim4_mp1;
  }

  template <>
  olim4_rhr & olim<olim4_rhr>() {
    return _olim._olim4_rhr;
  }

  template <>
  olim8_mp0 & olim<olim8_mp0>() {
    return _olim._olim8_mp0;
  }

  template <>
  olim8_mp1 & olim<olim8_mp1>() {
    return _olim._olim8_mp1;
  }

  template <>
  olim8_rhr & olim<olim8_rhr>() {
    return _olim._olim8_rhr;
  }

  template <>
  olim6_mp0 & olim<olim6_mp0>() {
    return _olim._olim6_mp0;
  }

  template <>
  olim6_mp1 & olim<olim6_mp1>() {
    return _olim._olim6_mp1;
  }

  template <>
  olim6_rhr & olim<olim6_rhr>() {
    return _olim._olim6_rhr;
  }

  template <>
  olim18_mp0 & olim<olim18_mp0>() {
    return _olim._olim18_mp0;
  }

  template <>
  olim18_mp1 & olim<olim18_mp1>() {
    return _olim._olim18_mp1;
  }

  template <>
  olim18_rhr & olim<olim18_rhr>() {
    return _olim._olim18_rhr;
  }

  template <>
  olim26_mp0 & olim<olim26_mp0>() {
    return _olim._olim26_mp0;
  }

  template <>
  olim26_mp1 & olim<olim26_mp1>() {
    return _olim._olim26_mp1;
  }

  template <>
  olim26_rhr & olim<olim26_rhr>() {
    return _olim._olim26_rhr;
  }

  template <>
  olim3d_mp0 & olim<olim3d_mp0>() {
    return _olim._olim3d_mp0;
  }

  template <>
  olim3d_mp1 & olim<olim3d_mp1>() {
    return _olim._olim3d_mp1;
  }

  template <>
  olim3d_rhr & olim<olim3d_rhr>() {
    return _olim._olim3d_rhr;
  }
};

template <class olim_t>
void construct_olim(olim_wrapper * w, olim_wrapper_params_s * p)
{
  new (&w->olim<olim_t>()) olim_t {{p->dims}, p->h, no_slow_t {}};
}

template <class olim_t>
void destruct_olim(olim_wrapper * w)
{
  w->olim<olim_t>().~olim_t();
}

status_e olim_wrapper_init(olim_wrapper ** w_ptr, olim_wrapper_params_s * p)
{
  olim_wrapper * w = (*w_ptr = new olim_wrapper {p->nb, p->F, {}});

  if (p->nb == OLIM4) {
    if (p->F == MP0) {
      construct_olim<olim4_mp0>(w, p);
    } else if (p->F == MP1) {
      construct_olim<olim4_mp1>(w, p);
    } else if (p->F == RHR) {
      construct_olim<olim4_rhr>(w, p);
    }
  } else if (p->nb == OLIM8) {
    if (p->F == MP0) {
      construct_olim<olim8_mp0>(w, p);
    } else if (p->F == MP1) {
      construct_olim<olim8_mp1>(w, p);
    } else if (p->F == RHR) {
      construct_olim<olim8_rhr>(w, p);
    }
  } else if (p->nb == OLIM6) {
    if (p->F == MP0) {
      construct_olim<olim6_mp0>(w, p);
    } else if (p->F == MP1) {
      construct_olim<olim6_mp1>(w, p);
    } else if (p->F == RHR) {
      construct_olim<olim6_rhr>(w, p);
    }
  } else if (p->nb == OLIM18) {
    if (p->F == MP0) {
      construct_olim<olim18_mp0>(w, p);
    } else if (p->F == MP1) {
      construct_olim<olim18_mp1>(w, p);
    } else if (p->F == RHR) {
      construct_olim<olim18_rhr>(w, p);
    }
  } else if (p->nb == OLIM26) {
    if (p->F == MP0) {
      construct_olim<olim26_mp0>(w, p);
    } else if (p->F == MP1) {
      construct_olim<olim26_mp1>(w, p);
    } else if (p->F == RHR) {
      construct_olim<olim26_rhr>(w, p);
    }
  } else if (p->nb == OLIM3D) {
    if (p->F == MP0) {
      construct_olim<olim3d_mp0>(w, p);
    } else if (p->F == MP1) {
      construct_olim<olim3d_mp1>(w, p);
    } else if (p->F == RHR) {
      construct_olim<olim3d_rhr>(w, p);
    }
  }

  return SUCCESS;
}

status_e olim_wrapper_deinit(olim_wrapper ** w_ptr)
{
  olim_wrapper * w = *w_ptr;

  if (w->nb == OLIM4) {
    if (w->F == MP0) {
      destruct_olim<olim4_mp0>(w);
    } else if (w->F == MP1) {
      destruct_olim<olim4_mp1>(w);
    } else if (w->F == RHR) {
      destruct_olim<olim4_rhr>(w);
    }
  } else if (w->nb == OLIM8) {
    if (w->F == MP0) {
      destruct_olim<olim8_mp0>(w);
    } else if (w->F == MP1) {
      destruct_olim<olim8_mp1>(w);
    } else if (w->F == RHR) {
      destruct_olim<olim8_rhr>(w);
    }
  } else if (w->nb == OLIM6) {
    if (w->F == MP0) {
      destruct_olim<olim6_mp0>(w);
    } else if (w->F == MP1) {
      destruct_olim<olim6_mp1>(w);
    } else if (w->F == RHR) {
      destruct_olim<olim6_rhr>(w);
    }
  } else if (w->nb == OLIM18) {
    if (w->F == MP0) {
      destruct_olim<olim18_mp0>(w);
    } else if (w->F == MP1) {
      destruct_olim<olim18_mp1>(w);
    } else if (w->F == RHR) {
      destruct_olim<olim18_rhr>(w);
    }
  } else if (w->nb == OLIM26) {
    if (w->F == MP0) {
      destruct_olim<olim26_mp0>(w);
    } else if (w->F == MP1) {
      destruct_olim<olim26_mp1>(w);
    } else if (w->F == RHR) {
      destruct_olim<olim26_rhr>(w);
    }
  } else if (w->nb == OLIM3D) {
    if (w->F == MP0) {
      destruct_olim<olim3d_mp0>(w);
    } else if (w->F == MP1) {
      destruct_olim<olim3d_mp1>(w);
    } else if (w->F == RHR) {
      destruct_olim<olim3d_rhr>(w);
    }
  }

  delete *w_ptr;

  return SUCCESS;
}

status_e olim_wrapper_run(olim_wrapper * w)
{
  if (w->nb == OLIM4) {
    if (w->F == MP0) {
      w->olim<olim4_mp0>().run();
    } else if (w->F == MP1) {
      w->olim<olim4_mp1>().run();
    } else if (w->F == RHR) {
      w->olim<olim4_rhr>().run();
    }
  } else if (w->nb == OLIM8) {
    if (w->F == MP0) {
      w->olim<olim8_mp0>().run();
    } else if (w->F == MP1) {
      w->olim<olim8_mp1>().run();
    } else if (w->F == RHR) {
      w->olim<olim8_rhr>().run();
    }
  } else if (w->nb == OLIM6) {
    if (w->F == MP0) {
      w->olim<olim6_mp0>().run();
    } else if (w->F == MP1) {
      w->olim<olim6_mp1>().run();
    } else if (w->F == RHR) {
      w->olim<olim6_rhr>().run();
    }
  } else if (w->nb == OLIM18) {
    if (w->F == MP0) {
      w->olim<olim18_mp0>().run();
    } else if (w->F == MP1) {
      w->olim<olim18_mp1>().run();
    } else if (w->F == RHR) {
      w->olim<olim18_rhr>().run();
    }
  } else if (w->nb == OLIM26) {
    if (w->F == MP0) {
      w->olim<olim26_mp0>().run();
    } else if (w->F == MP1) {
      w->olim<olim26_mp1>().run();
    } else if (w->F == RHR) {
      w->olim<olim26_rhr>().run();
    }
  } else if (w->nb == OLIM3D) {
    if (w->F == MP0) {
      w->olim<olim3d_mp0>().run();
    } else if (w->F == MP1) {
      w->olim<olim3d_mp1>().run();
    } else if (w->F == RHR) {
      w->olim<olim3d_rhr>().run();
    }
  }

  return SUCCESS;
}

status_e olim_wrapper_add_src(olim_wrapper * w, int * inds, double U)
{
  if (w->nb == OLIM4) {
    if (w->F == MP0) {
      w->olim<olim4_mp0>().add_src(inds, U);
    } else if (w->F == MP1) {
      w->olim<olim4_mp1>().add_src(inds, U);
    } else if (w->F == RHR) {
      w->olim<olim4_rhr>().add_src(inds, U);
    }
  } else if (w->nb == OLIM8) {
    if (w->F == MP0) {
      w->olim<olim8_mp0>().add_src(inds, U);
    } else if (w->F == MP1) {
      w->olim<olim8_mp1>().add_src(inds, U);
    } else if (w->F == RHR) {
      w->olim<olim8_rhr>().add_src(inds, U);
    }
  } else if (w->nb == OLIM6) {
    if (w->F == MP0) {
      w->olim<olim6_mp0>().add_src(inds, U);
    } else if (w->F == MP1) {
      w->olim<olim6_mp1>().add_src(inds, U);
    } else if (w->F == RHR) {
      w->olim<olim6_rhr>().add_src(inds, U);
    }
  } else if (w->nb == OLIM18) {
    if (w->F == MP0) {
      w->olim<olim18_mp0>().add_src(inds, U);
    } else if (w->F == MP1) {
      w->olim<olim18_mp1>().add_src(inds, U);
    } else if (w->F == RHR) {
      w->olim<olim18_rhr>().add_src(inds, U);
    }
  } else if (w->nb == OLIM26) {
    if (w->F == MP0) {
      w->olim<olim26_mp0>().add_src(inds, U);
    } else if (w->F == MP1) {
      w->olim<olim26_mp1>().add_src(inds, U);
    } else if (w->F == RHR) {
      w->olim<olim26_rhr>().add_src(inds, U);
    }
  } else if (w->nb == OLIM3D) {
    if (w->F == MP0) {
      w->olim<olim3d_mp0>().add_src(inds, U);
    } else if (w->F == MP1) {
      w->olim<olim3d_mp1>().add_src(inds, U);
    } else if (w->F == RHR) {
      w->olim<olim3d_rhr>().add_src(inds, U);
    }
  }

  return SUCCESS;
}

status_e olim_wrapper_add_bd(olim_wrapper * w, int * inds)
{
  if (w->nb == OLIM4) {
    if (w->F == MP0) {
      w->olim<olim4_mp0>().add_bd(inds);
    } else if (w->F == MP1) {
      w->olim<olim4_mp1>().add_bd(inds);
    } else if (w->F == RHR) {
      w->olim<olim4_rhr>().add_bd(inds);
    }
  } else if (w->nb == OLIM8) {
    if (w->F == MP0) {
      w->olim<olim8_mp0>().add_bd(inds);
    } else if (w->F == MP1) {
      w->olim<olim8_mp1>().add_bd(inds);
    } else if (w->F == RHR) {
      w->olim<olim8_rhr>().add_bd(inds);
    }
  } else if (w->nb == OLIM6) {
    if (w->F == MP0) {
      w->olim<olim6_mp0>().add_bd(inds);
    } else if (w->F == MP1) {
      w->olim<olim6_mp1>().add_bd(inds);
    } else if (w->F == RHR) {
      w->olim<olim6_rhr>().add_bd(inds);
    }
  } else if (w->nb == OLIM18) {
    if (w->F == MP0) {
      w->olim<olim18_mp0>().add_bd(inds);
    } else if (w->F == MP1) {
      w->olim<olim18_mp1>().add_bd(inds);
    } else if (w->F == RHR) {
      w->olim<olim18_rhr>().add_bd(inds);
    }
  } else if (w->nb == OLIM26) {
    if (w->F == MP0) {
      w->olim<olim26_mp0>().add_bd(inds);
    } else if (w->F == MP1) {
      w->olim<olim26_mp1>().add_bd(inds);
    } else if (w->F == RHR) {
      w->olim<olim26_rhr>().add_bd(inds);
    }
  } else if (w->nb == OLIM3D) {
    if (w->F == MP0) {
      w->olim<olim3d_mp0>().add_bd(inds);
    } else if (w->F == MP1) {
      w->olim<olim3d_mp1>().add_bd(inds);
    } else if (w->F == RHR) {
      w->olim<olim3d_rhr>().add_bd(inds);
    }
  }

  return SUCCESS;
}

status_e olim_wrapper_get_U_ptr(olim_wrapper * w, double ** U_ptr)
{
  if (w->nb == OLIM4) {
    if (w->F == MP0) {
      *U_ptr = w->olim<olim4_mp0>().get_U_ptr();
    } else if (w->F == MP1) {
      *U_ptr = w->olim<olim4_mp1>().get_U_ptr();
    } else if (w->F == RHR) {
      *U_ptr = w->olim<olim4_rhr>().get_U_ptr();
    }
  } else if (w->nb == OLIM8) {
    if (w->F == MP0) {
      *U_ptr = w->olim<olim8_mp0>().get_U_ptr();
    } else if (w->F == MP1) {
      *U_ptr = w->olim<olim8_mp1>().get_U_ptr();
    } else if (w->F == RHR) {
      *U_ptr = w->olim<olim8_rhr>().get_U_ptr();
    }
  } else if (w->nb == OLIM6) {
    if (w->F == MP0) {
      *U_ptr = w->olim<olim6_mp0>().get_U_ptr();
    } else if (w->F == MP1) {
      *U_ptr = w->olim<olim6_mp1>().get_U_ptr();
    } else if (w->F == RHR) {
      *U_ptr = w->olim<olim6_rhr>().get_U_ptr();
    }
  } else if (w->nb == OLIM18) {
    if (w->F == MP0) {
      *U_ptr = w->olim<olim18_mp0>().get_U_ptr();
    } else if (w->F == MP1) {
      *U_ptr = w->olim<olim18_mp1>().get_U_ptr();
    } else if (w->F == RHR) {
      *U_ptr = w->olim<olim18_rhr>().get_U_ptr();
    }
  } else if (w->nb == OLIM26) {
    if (w->F == MP0) {
      *U_ptr = w->olim<olim26_mp0>().get_U_ptr();
    } else if (w->F == MP1) {
      *U_ptr = w->olim<olim26_mp1>().get_U_ptr();
    } else if (w->F == RHR) {
      *U_ptr = w->olim<olim26_rhr>().get_U_ptr();
    }
  } else if (w->nb == OLIM3D) {
    if (w->F == MP0) {
      *U_ptr = w->olim<olim3d_mp0>().get_U_ptr();
    } else if (w->F == MP1) {
      *U_ptr = w->olim<olim3d_mp1>().get_U_ptr();
    } else if (w->F == RHR) {
      *U_ptr = w->olim<olim3d_rhr>().get_U_ptr();
    }
  }

  return SUCCESS;
}

status_e olim_wrapper_get_s_ptr(olim_wrapper * w, double ** s_ptr)
{
  if (w->nb == OLIM4) {
    if (w->F == MP0) {
      *s_ptr = w->olim<olim4_mp0>().get_s_ptr();
    } else if (w->F == MP1) {
      *s_ptr = w->olim<olim4_mp1>().get_s_ptr();
    } else if (w->F == RHR) {
      *s_ptr = w->olim<olim4_rhr>().get_s_ptr();
    }
  } else if (w->nb == OLIM8) {
    if (w->F == MP0) {
      *s_ptr = w->olim<olim8_mp0>().get_s_ptr();
    } else if (w->F == MP1) {
      *s_ptr = w->olim<olim8_mp1>().get_s_ptr();
    } else if (w->F == RHR) {
      *s_ptr = w->olim<olim8_rhr>().get_s_ptr();
    }
  } else if (w->nb == OLIM6) {
    if (w->F == MP0) {
      *s_ptr = w->olim<olim6_mp0>().get_s_ptr();
    } else if (w->F == MP1) {
      *s_ptr = w->olim<olim6_mp1>().get_s_ptr();
    } else if (w->F == RHR) {
      *s_ptr = w->olim<olim6_rhr>().get_s_ptr();
    }
  } else if (w->nb == OLIM18) {
    if (w->F == MP0) {
      *s_ptr = w->olim<olim18_mp0>().get_s_ptr();
    } else if (w->F == MP1) {
      *s_ptr = w->olim<olim18_mp1>().get_s_ptr();
    } else if (w->F == RHR) {
      *s_ptr = w->olim<olim18_rhr>().get_s_ptr();
    }
  } else if (w->nb == OLIM26) {
    if (w->F == MP0) {
      *s_ptr = w->olim<olim26_mp0>().get_s_ptr();
    } else if (w->F == MP1) {
      *s_ptr = w->olim<olim26_mp1>().get_s_ptr();
    } else if (w->F == RHR) {
      *s_ptr = w->olim<olim26_rhr>().get_s_ptr();
    }
  } else if (w->nb == OLIM3D) {
    if (w->F == MP0) {
      *s_ptr = w->olim<olim3d_mp0>().get_s_ptr();
    } else if (w->F == MP1) {
      *s_ptr = w->olim<olim3d_mp1>().get_s_ptr();
    } else if (w->F == RHR) {
      *s_ptr = w->olim<olim3d_rhr>().get_s_ptr();
    }
  }

  return SUCCESS;
}

status_e olim_wrapper_get_state_ptr(olim_wrapper * w, char ** state_ptr)
{
  if (w->nb == OLIM4) {
    if (w->F == MP0) {
      *state_ptr = w->olim<olim4_mp0>().get_state_ptr();
    } else if (w->F == MP1) {
      *state_ptr = w->olim<olim4_mp1>().get_state_ptr();
    } else if (w->F == RHR) {
      *state_ptr = w->olim<olim4_rhr>().get_state_ptr();
    }
  } else if (w->nb == OLIM8) {
    if (w->F == MP0) {
      *state_ptr = w->olim<olim8_mp0>().get_state_ptr();
    } else if (w->F == MP1) {
      *state_ptr = w->olim<olim8_mp1>().get_state_ptr();
    } else if (w->F == RHR) {
      *state_ptr = w->olim<olim8_rhr>().get_state_ptr();
    }
  } else if (w->nb == OLIM6) {
    if (w->F == MP0) {
      *state_ptr = w->olim<olim6_mp0>().get_state_ptr();
    } else if (w->F == MP1) {
      *state_ptr = w->olim<olim6_mp1>().get_state_ptr();
    } else if (w->F == RHR) {
      *state_ptr = w->olim<olim6_rhr>().get_state_ptr();
    }
  } else if (w->nb == OLIM18) {
    if (w->F == MP0) {
      *state_ptr = w->olim<olim18_mp0>().get_state_ptr();
    } else if (w->F == MP1) {
      *state_ptr = w->olim<olim18_mp1>().get_state_ptr();
    } else if (w->F == RHR) {
      *state_ptr = w->olim<olim18_rhr>().get_state_ptr();
    }
  } else if (w->nb == OLIM26) {
    if (w->F == MP0) {
      *state_ptr = w->olim<olim26_mp0>().get_state_ptr();
    } else if (w->F == MP1) {
      *state_ptr = w->olim<olim26_mp1>().get_state_ptr();
    } else if (w->F == RHR) {
      *state_ptr = w->olim<olim26_rhr>().get_state_ptr();
    }
  } else if (w->nb == OLIM3D) {
    if (w->F == MP0) {
      *state_ptr = w->olim<olim3d_mp0>().get_state_ptr();
    } else if (w->F == MP1) {
      *state_ptr = w->olim<olim3d_mp1>().get_state_ptr();
    } else if (w->F == RHR) {
      *state_ptr = w->olim<olim3d_rhr>().get_state_ptr();
    }
  }

  return SUCCESS;
}
