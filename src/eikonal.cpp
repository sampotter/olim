#include "eikonal.h"

#include <variant>

#include "eikonal/fmm.hpp"
#include "eikonal/olim.hpp"
#include "eikonal/olim3d.hpp"

namespace eikonal {

template <ordering ord = ordering::ROW_MAJOR>
using olim_v = std::variant<
  fmm2<ord>,
  fmm3<ord>,
  olim4_mp0<ord>,
  olim4_mp1<ord>,
  olim4_rhr<ord>,
  olim8_mp0<ord>,
  olim8_mp1<ord>,
  olim8_rhr<ord>,
  olim6_mp0<ord>,
  olim6_mp1<ord>,
  olim6_rhr<ord>,
  olim18_mp0<ord>,
  olim18_mp1<ord>,
  olim18_rhr<ord>,
  olim26_mp0<ord>,
  olim26_mp1<ord>,
  olim26_rhr<ord>,
  olim3d_mp0<ord>,
  olim3d_mp1<ord>,
  olim3d_rhr<ord>
>;

}

struct eikonal_wkspc
{
  static constexpr ordering ord = ordering::ROW_MAJOR;

  eikonal_wkspc(eikonal::olim_v<ord> olim): _olim {olim} {}

  ~eikonal_wkspc() {
    // TODO: ...
  }

  eikonal::olim_v<ord> _olim;
};

template <class olim_t>
void construct_olim_variant(eikonal_wkspc ** w_ptr, eikonal_params_s * p) {
  *w_ptr = new eikonal_wkspc {olim_t {{p->dims}, p->h, no_slow_t {}}};
}

status_e eikonal_init(eikonal_wkspc ** w_ptr, eikonal_params_s * p)
{
  constexpr ordering ord = eikonal_wkspc::ord;

  if (p->nb == OLIM4) {
    if (p->F == MP0) {
      construct_olim_variant<eikonal::olim4_mp0<ord>>(w_ptr, p);
    } else if (p->F == MP1) {
      construct_olim_variant<eikonal::olim4_mp1<ord>>(w_ptr, p);
    } else if (p->F == RHR) {
      construct_olim_variant<eikonal::olim4_rhr<ord>>(w_ptr, p);
    }
  } else if (p->nb == OLIM8) {
    if (p->F == MP0) {
      construct_olim_variant<eikonal::olim8_mp0<ord>>(w_ptr, p);
    } else if (p->F == MP1) {
      construct_olim_variant<eikonal::olim8_mp1<ord>>(w_ptr, p);
    } else if (p->F == RHR) {
      construct_olim_variant<eikonal::olim8_rhr<ord>>(w_ptr, p);
    }
  } else if (p->nb == OLIM6) {
    if (p->F == MP0) {
      construct_olim_variant<eikonal::olim6_mp0<ord>>(w_ptr, p);
    } else if (p->F == MP1) {
      construct_olim_variant<eikonal::olim6_mp1<ord>>(w_ptr, p);
    } else if (p->F == RHR) {
      construct_olim_variant<eikonal::olim6_rhr<ord>>(w_ptr, p);
    }
  } else if (p->nb == OLIM18) {
    if (p->F == MP0) {
      construct_olim_variant<eikonal::olim18_mp0<ord>>(w_ptr, p);
    } else if (p->F == MP1) {
      construct_olim_variant<eikonal::olim18_mp1<ord>>(w_ptr, p);
    } else if (p->F == RHR) {
      construct_olim_variant<eikonal::olim18_rhr<ord>>(w_ptr, p);
    }
  } else if (p->nb == OLIM26) {
    if (p->F == MP0) {
      construct_olim_variant<eikonal::olim26_mp0<ord>>(w_ptr, p);
    } else if (p->F == MP1) {
      construct_olim_variant<eikonal::olim26_mp1<ord>>(w_ptr, p);
    } else if (p->F == RHR) {
      construct_olim_variant<eikonal::olim26_rhr<ord>>(w_ptr, p);
    }
  } else if (p->nb == OLIM3D) {
    if (p->F == MP0) {
      construct_olim_variant<eikonal::olim3d_mp0<ord>>(w_ptr, p);
    } else if (p->F == MP1) {
      construct_olim_variant<eikonal::olim3d_mp1<ord>>(w_ptr, p);
    } else if (p->F == RHR) {
      construct_olim_variant<eikonal::olim3d_rhr<ord>>(w_ptr, p);
    }
  } else if (p->nb == FMM2) {
    if (p->F == RHR) {
      construct_olim_variant<eikonal::fmm<2, ord>>(w_ptr, p);
    } else {
      throw std::runtime_error("FMM2 requires quad == RHR");
    }
  } else if (p->nb == FMM3) {
    if (p->F == RHR) {
      construct_olim_variant<eikonal::fmm<3, ord>>(w_ptr, p);
    } else {
      throw std::runtime_error("FMM3 requires quad == RHR");
    }
  }

  return SUCCESS;
}

status_e eikonal_deinit(eikonal_wkspc ** w_ptr)
{
  delete *w_ptr;

  *w_ptr = nullptr;

  return SUCCESS;
}

status_e eikonal_solve(eikonal_wkspc * w)
{
  std::visit([] (auto && o) { o.solve(); }, w->_olim);

  return SUCCESS;
}

status_e eikonal_step(eikonal_wkspc * w, int *lin)
{
  std::visit(
    [&] (auto && o) { *lin = o.step(); },
    w->_olim
  );

  return SUCCESS;
}

status_e eikonal_adjust(eikonal_wkspc * w, int * inds, double U)
{
  std::visit(
    [&] (auto && o) { o.adjust(inds, U); },
    w->_olim
  );

  return SUCCESS;
}

status_e eikonal_add_src(eikonal_wkspc * w, int * inds, double U)
{
  std::visit(
    /**
     * XXX: the next line doesn't compile with gcc... need to check if
     * the following uncommented line *does* compile with clang
     */
    // [&] (auto && o) { o.base_marcher::add_src(inds, U); },
    [&] (auto && o) { o.add_src(inds, U); },
    w->_olim
  );

  return SUCCESS;
}

status_e
eikonal_factor(eikonal_wkspc * w, int * inds, fac_src_s * fs)
{
  std::visit(
    [&] (auto && o) { o.factor(inds, fs); },
    w->_olim
  );

  return SUCCESS;
}

status_e eikonal_add_bd(eikonal_wkspc * w, int * inds)
{
  std::visit(
    [&] (auto && o) { o.add_bd(inds); },
    w->_olim
  );

  return SUCCESS;
}

status_e eikonal_peek(eikonal_wkspc * w, double * value, int * lin, bool * empty)
{
  std::visit(
    [&] (auto && o) { *empty = o.peek(value, lin); },
    w->_olim
  );

  return SUCCESS;
}

status_e eikonal_get_U_ptr(eikonal_wkspc * w, double ** U_ptr)
{
  std::visit(
    [&] (auto && o) { *U_ptr = o.get_U_ptr(); },
    w->_olim
  );

  return SUCCESS;
}

status_e eikonal_get_s_ptr(eikonal_wkspc * w, double ** s_ptr)
{
  std::visit(
    [&] (auto && o) { *s_ptr = o.get_s_ptr(); },
    w->_olim
  );

  return SUCCESS;
}

status_e eikonal_get_state_ptr(eikonal_wkspc * w, char ** state_ptr)
{
  std::visit(
    [&] (auto && o) { *state_ptr = o.get_state_ptr(); },
    w->_olim
  );


  return SUCCESS;
}
