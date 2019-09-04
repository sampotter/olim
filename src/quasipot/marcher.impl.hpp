#pragma once

#include "../mat.hpp"

template <class derived, int n, ordering ord>
quasipot::marcher<derived, n, ord>::marcher(ivec dims, double h, vfield b, int K):
  base_marcher<marcher<derived, n, ord>, n> {dims},
  _h {h},
  _b {b},
  _K {K},
  _update_box_dims {(2*K + 1)*ivec::one()},
  _center {K*ivec::one()},
  _valid_front {(2*K + 3)*ivec::one()},
  _valid_front_linear_indices {new int[pow(2*K + 1, n)]}
{
  assert(_K > 0);
}

template <class derived, int n, ordering ord>
quasipot::marcher<derived, n, ord>::~marcher()
{
  delete[] _valid_front_linear_indices;
}

template <class derived, int n, ordering ord>
void
quasipot::marcher<derived, n, ord>::visit_neighbors(int lin_center)
{
  set_valid_front(lin_center);

  auto const update = [&] (int lin, state s) {
    auto U = inf<double>;
    static_cast<derived *>(this)->update_impl(lin, U);
    if (U < this->_U[lin]) {
      this->_U[lin] = U;
      if (s == state::trial) {
        this->_heap.update(lin);
      }
    }
  };

  for (int i = 0; i < get_num_nb(); ++i) {
    int lin = lin_center + this->_linear_offset[i];
    auto s = this->_state[lin];
    if (s == state::trial) {
      set_valid_front_linear_indices(lin);
      update(lin, s);
    }
  }
}

template <class derived, int n, ordering ord>
bool
quasipot::marcher<derived, n, ord>::is_valid_front(ivec inds) const
{
  for (auto offset: range<n, ord> {3*ivec::one()}) {
    ivec nb_inds = inds + offset - ivec::one();
    state s = this->get_state(nb_inds);
    if (this->in_bounds(nb_inds) && (s == state::trial || s == state::far)) {
      return true;
    }
  }
  return false;
}

template <class derived, int n, ordering ord>
void
quasipot::marcher<derived, n, ord>::set_valid_front(int lin)
{
  auto inds = this->to_vector_index(lin);

  _valid_front.fill(false);
  for (auto offset: range<n, ord> {_valid_front.get_dims()}) {
    _valid_front(offset) = is_valid_front(inds - _center + offset);
  }
}

template <class derived, int n, ordering ord>
void
quasipot::marcher<derived, n, ord>::set_valid_front_linear_indices(int lin)
{
  ivec inds = this->to_vector_index(lin);

  ivec box_offset = inds - _newly_valid;
  ivec lo = ivec::one() + box_offset;
  ivec hi = lo + _update_box_dims;
  auto valid_front_subgrid = _valid_front.subgrid(lo, hi);

  int size = 0;
  for (auto offset: range<n, ord> {_update_box_dims}) {
    if (valid_front_subgrid(offset)) {
      _valid_front_linear_indices[size++] =
        this->to_linear_index(inds - _center + box_offset + offset);
    }
  }
}
