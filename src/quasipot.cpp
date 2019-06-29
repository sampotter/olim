#include "quasipot.h"

#include <variant>

#include "quasipot/olim.hpp"

using quasipot_v = std::variant<quasipot::olim, quasipot::olim3d>;

struct quasipot
{
  quasipot_v * impl;
};

status_e quasipot_init(quasipot_s **handle, quasipot_params_s *p)
{
  *handle = new quasipot;

  (*handle)->impl = p->ndims == 2 ?
    new quasipot_v(
      std::in_place_type_t<quasipot::olim>,
      {p->dims[0], p->dims[1]},
      p->h,
      (quasipot::olim::vfield) p->b,
      p->K):
    new quasipot_v(
      std::in_place_type_t<quasipot::olim3d>,
      {p->dims[0], p->dims[1], p->dims[2]},
      p->h,
      (quasipot::olim::vfield) p->b,
      p->K);

  return SUCCESS;
}

status_e quasipot_deinit(quasipot_s **handle)
{
  if (*handle == nullptr) {
    return SUCCESS;
  }

  delete ((*handle)->impl);
  delete *handle;

  return SUCCESS;
}

status_e quasipot_solve(quasipot_s *q)
{
  std::visit(
    [] (auto & _) {
      _.solve();
    },
    *q->impl);
}

status_e quasipot_step(quasipot_s *q, int *lin)
{
  std::visit(
    [] (auto & _) {
      *lin = _.step();
    },
    *q->impl);
}

status_e quasipot_peek(quasipot_s *q, double *value, int *lin, bool *empty)
{
  std::visit(
    [] (auto & _) {
      *empty = _.peek(value, lin);
    },
    *q->impl);
}
