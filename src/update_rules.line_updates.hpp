#ifndef __UPDATE_RULES_LINE_UPDATES_HPP__
#define __UPDATE_RULES_LINE_UPDATES_HPP__

namespace update_rules {
  struct rhr_line_updates {
    template <int D>
    double line(double u0, double s, double s0, double h) const;
  };

  struct mp_line_updates {
    template <int D>
    double line(double u0, double s, double s0, double h) const;
  };
}

#include "update_rules.line_updates.impl.hpp"

#endif // __UPDATE_RULES_LINE_UPDATES_HPP__
