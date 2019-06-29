#pragma once

#include "../grid.hpp"

namespace quasipot {

template<class base, int n, ordering ord>
struct marcher: base_marcher<marcher<base, n, ord>, n>
{
  using fvec = vec<double, n>;
  using ivec = vec<int, n>;
  using vfield = fvec(*)(fvec);

  static constexpr int ndim = n; // TODO: -> get_n()

  static constexpr int get_nnb() { return detail::max_num_nb(); }
  static constexpr ordering get_ord() { return ord; }

  marcher(ivec dims, double h, vfield b, int K);
  ~marcher();

  void add_src(int * inds);
  void add_src(ivec inds);

OLIM_PROTECTED:

  inline double get_h() const { return _h; }
  void visit_neighbors(int lin);

  double _h;
  vfield _b;
  int _K;

OLIM_PRIVATE:

  bool is_valid_front(ivec inds) const;
  void set_valid_front(ivec inds);
  void set_valid_front_linear_indices(ivec inds);

  inline fvec get_x(ivec inds1, ivec inds0) const {
    return _h*static_cast<fvec>(inds1 - inds0);
  }

  // TODO: we could try using an LRU cache to reduce the amount we have
  // to update `_valid_front`...
  //
  // We would probably need to introduce a "tristate", but we could
  // have the grid be "unset/true/false". Then, in the LRU cache, we
  // keep track of some number of recently used grids, and initialize
  // a new grid using any of the grids that intersect the current grid

  ivec _center;
  ivec _newly_valid;

  /**
   * Binary mask indicating which nodes in the (2K+3)^n box centered
   * around the new valid node are in the "valid front" (i.e., are
   * valid and on the boundary of the set of valid nodes).
   *
   * This gets reset each time a node comes off the heap.
   *
   * The box is (2K+3)^n because we use subgrid to extract the relevant
   * "valid front" nodes when we update each trial node that neighbors
   * the new valid node.
   */
  grid<bool, n, ord> _valid_front;

  /**
   * An array containing a list of linear indices to the "valid front" nodes
   * in the (2K+1)^n box surrounding a trial node that is being updated.
   *
   * This gets updated before we start updating one of the trial nodes
   * surrounding the new valid node and is a subset of the indices in
   * `_valid_front`.
   */
  int * _valid_front_linear_indices;
  int _valid_front_size; // current number of indices
};

#include "marcher.impl.hpp"
