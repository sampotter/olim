#ifndef __UPDATE_RULES_LINE_UPDATES_HPP__
#define __UPDATE_RULES_LINE_UPDATES_HPP__

namespace update_rules {
  struct line_updates {
    /**
     * This line update is for use with "rectangular" OLIMs (RHR and
     * MP0).
     */
    template <int N>
    double line(double u0, double s, double h) const;
    
    /**
     * This line update is for the MP1 OLIMs.
     *
     * TODO: the name clash between this line update function and the
     * preceding function is probably not a big deal, but if we ever
     * start looking at methods which use a higher-order line update,
     * then we'll need to come up with a better naming or template
     * specialization convention.
     */
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
