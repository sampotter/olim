#ifndef __UPDATE_RULES_LINE_UPDATES_HPP__
#define __UPDATE_RULES_LINE_UPDATES_HPP__

namespace update_rules {
  struct line_updates {
    template <int N>
    double line(double u0, double s, double h) const;
    
    template <int N>
    double line(double u0, double s, double s0, double h) const;
  };
}

#include "update_rules.line_updates.impl.hpp"

#endif // __UPDATE_RULES_LINE_UPDATES_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
