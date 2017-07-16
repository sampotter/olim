#ifndef __SOLIM4_NODE_HPP__
#define __SOLIM4_NODE_HPP__

#include "smart_node.hpp"

struct solim4_node: public smart_node {
  enum class update_type {ONE_POINT, TWO_POINT};

  update_type get_update_type() const { return _update_type; }
  void set_update_type(update_type t) { _update_type = t; }
private:
  update_type _update_type;
};

#endif // __SOLIM4_NODE_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
