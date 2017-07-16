#include "abstract_marcher.hpp"

abstract_marcher::abstract_marcher(double h, int S_size):
  _h {h},
  _S_cache {new double[S_size]},
  _should_free_S_cache {true}
{
  for (int i = 0; i < S_size; ++i) {
    _S_cache[i] = -1;
  }
}

abstract_marcher::abstract_marcher(double h, double * S_cache):
  _h {h},
  _S_cache {S_cache},
  _should_free_S_cache {false}
{}

abstract_marcher::~abstract_marcher()
{
  if (_should_free_S_cache) {
    delete[] _S_cache;
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:

