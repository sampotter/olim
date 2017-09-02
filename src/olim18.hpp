#ifndef __OLIM18_HPP__
#define __OLIM18_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"

template <class Node, class Updates>
struct olim18: public marcher_3d<Node>, public Updates {
  using marcher_3d<Node>::marcher_3d;
protected:
  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb);
  static int di[6];
  static int dj[6];
  static int dk[6];
private:
  virtual void stage_neighbors_impl(abstract_node * n);
  virtual void update_impl(int i, int j, int k, double & T);
};

#endif // __OLIM18_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
