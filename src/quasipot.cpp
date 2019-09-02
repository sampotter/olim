#include "quasipot.h"

#include <variant>

#include "quasipot/olim.hpp"

namespace quasipot {

template <ordering ord>
using olim_v = std::variant<
  olim<ord>
>;

}

struct quasipot_wkspc
{
  static constexpr ordering ord = ordering::ROW_MAJOR;

  quasipot_wkspc(quasipot::olim_v<ord> olim): _olim {olim} {}

  quasipot::olim_v<ord> _olim;
};

status_e quasipot_init(quasipot_wkspc_s **w_ptr, quasipot_params_s *p)
{
  constexpr ordering ord = quasipot_wkspc::ord;

  *w_ptr = new quasipot_wkspc {
    quasipot::olim<ord> {{p->dims}, p->h, p->b, p->K}
  };

  return SUCCESS;
}

status_e quasipot_deinit(quasipot_wkspc_s **w_ptr)
{
  if (*w_ptr == nullptr) {
    return SUCCESS;
  }

  delete *w_ptr;
  *w_ptr = nullptr;

  return SUCCESS;
}

status_e quasipot_solve(quasipot_wkspc_s *w)
{
  std::visit(
    [] (auto & _) { _.solve(); },
    w->_olim
  );

  return SUCCESS;
}

status_e quasipot_step(quasipot_wkspc_s *w, int *lin)
{
  std::visit(
    [&] (auto & _) { *lin = _.step(); },
    w->_olim
  );

  return SUCCESS;
}

status_e quasipot_peek(quasipot_wkspc_s *w, double *value, int *lin,
                       bool *empty)
{
  std::visit(
    [&] (auto & _) { *empty = _.peek(value, lin); },
    w->_olim
  );

  return SUCCESS;
}

status_e quasipot_adjust(quasipot_wkspc_s *w, int *inds, double U)
{
  std::visit(
    [&] (auto & _) { _.adjust(inds, U); },
    w->_olim
  );

  return SUCCESS;
}

status_e quasipot_add_src(quasipot_wkspc_s *w, int *inds, double U)
{
  std::visit(
    [&] (auto & _) { _.add_src(inds, U); },
    w->_olim
  );

  return SUCCESS;
}

status_e quasipot_add_bd(quasipot_wkspc_s *w, int *inds)
{
  std::visit(
    [&] (auto & _) { _.add_bd(inds); },
    w->_olim
  );

  return SUCCESS;
}

status_e quasipot_add_free(quasipot_wkspc_s *w, int *inds)
{
  std::visit(
    [&] (auto & _) { _.add_free(inds); },
    w->_olim
  );

  return SUCCESS;
}

status_e quasipot_get_U_ptr(quasipot_wkspc_s *w, double **U_ptr)
{
  std::visit(
    [&] (auto & _) { *U_ptr = _.get_U_ptr(); },
    w->_olim
  );

  return SUCCESS;
}

status_e quasipot_get_state_ptr(quasipot_wkspc_s *w, char **state_ptr)
{
  std::visit(
    [&] (auto & _) { *state_ptr = _.get_state_ptr(); },
    w->_olim
  );

  return SUCCESS;
}
