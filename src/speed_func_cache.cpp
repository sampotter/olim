#include "speed_func_cache.hpp"

speed_func_cache::speed_func_cache(int size, double * cache_ptr):
  _should_free {cache_ptr == nullptr},
  _cache {_should_free ? new double[size] : cache_ptr}
{
  // If the user didn't pass a pointer to an S_cache, we need to
  // initialize the one we've just created.
  if (_should_free) {
    for (int i = 0; i < size; ++i) {
      _cache[i] = -1;
    }
  }
}

speed_func_cache::~speed_func_cache()
{
  if (_should_free) {
    delete[] _cache;
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
