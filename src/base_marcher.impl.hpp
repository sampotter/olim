#pragma once

template <class base, int n>
void base_marcher<base, n>::solve()
{
  while (!_heap.empty()) {
    (void) step_impl();
  }
}

template <class base, int n>
int base_marcher<base, n>::step()
{
  return _heap.empty() ?
    -1 :
    static_cast<base *>(this)->to_external_linear_index(step_impl());
}

template <class base, int n>
int base_marcher<base, n>::step_impl()
{
  int lin = _heap.front();
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(_state[lin] == state::trial);
#endif
  _heap.pop_front();
  _state[lin] = state::valid;
  stage_neighbors(lin);
  static_cast<base *>(this)->visit_neighbors(lin);
  return lin;
}

template <class base, int n>
void base_marcher<base, n>::add_src(int const * inds, double U) {
  add_src(ivec {inds}, U);
}

template <class base, int n>
void base_marcher<base, n>::add_src(ivec inds, double U) {
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  inds += ivec::one();
  int lin = to_linear_index(inds);
  _U[lin] = U;
  _state[lin] = state::trial;
  _heap.insert(lin);
}

template <class base, int n>
void
base_marcher<base, n>::add_srcs(ivec const * inds, double const * U, int num) {
  for (int i = 0; i < num; ++i) {
    add_src(inds[i], U[i]);
  }
}

template <class base, int n>
void
base_marcher<base, n>::add_bd(int const * inds)
{
  add_bd(ivec {inds});
}

template <class base, int n>
void
base_marcher<base, n>::add_bd(ivec inds)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  inds += ivec::one();
  int lin = to_linear_index(inds);
  _U[lin] = inf<double>;
  _state[lin] = state::boundary;
}

template <class base, int n>
bool base_marcher<base, n>::peek(double * U, int * lin) const {
  bool const is_empty = _heap.empty();
  if (!is_empty) {
    int lin_ = _heap.front();
    if (U != nullptr) {
      *U = _U[lin_];
    }
    if (lin != nullptr) {
      *lin = static_cast<base const *>(this)->to_external_linear_index(lin_);
    }
  }
  return is_empty;
}

template <class base, int n>
void base_marcher<base, n>::adjust(int const * inds, double U) {
  adjust(ivec {inds}, U);
}

template <class base, int n>
void base_marcher<base, n>::adjust(ivec inds, double U) {
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(this->in_bounds(inds));
#endif
  inds += ivec::one();
  int lin = static_cast<base *>(this)->to_linear_index(inds);
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(U <= this->_U[lin]);
  assert(this->_state[lin] == state::trial);
#endif
  this->_U[lin] = U;
  this->_heap.update(lin);
}

template <class base, int n>
double base_marcher<base, n>::get_U(ivec inds) const {
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(this->_U != nullptr);
  assert(in_bounds(inds));
#endif
  return this->_U[
    static_cast<base const *>(this)->to_linear_index(inds + ivec::one())];
}

template <class base, int n>
state base_marcher<base, n>::get_state(ivec inds) const {
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(this->_state != nullptr);
  assert(in_bounds(inds));
#endif
  return this->_state[to_linear_index(inds + ivec::one())];
}

template <class base, int n>
void base_marcher<base, n>::stage_neighbors(int lin_center) {
  // Traverse the neighborhood of the node at `lin_center`, set all
  // far nodes to trial, and insert them into the heap.
  for (int i = 0; i < base::get_num_nb(); ++i) {
    int lin = lin_center + this->_linear_offset[i];
    if (this->_state[lin] == state::far) {
      this->_state[lin] = state::trial;
      this->_heap.insert(lin);
    }
  }
}

template <class base, int n>
bool base_marcher<base, n>::in_bounds(ivec inds) const {
  return uvec {inds} < uvec {this->_dims - 2*ivec::one()};
}
