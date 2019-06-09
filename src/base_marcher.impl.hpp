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
  static_cast<base *>(this)->visit_neighbors(lin);
  return lin;
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
