#include "marcher.hpp"

#include <cassert>

marcher::marcher(int height, int width, double h, speed_func S,
                 double x0, double y0):
  abstract_marcher {h, width*height},
  _S {S},
  _x0 {x0},
  _y0 {y0},
  _height {height},
  _width {width}
{}

marcher::marcher(int height, int width, double h, double * S_cache):
  abstract_marcher {h, width*height, S_cache},
  _height {height},
  _width {width}
{}

double marcher::S(int i, int j) {
  assert(in_bounds(i, j));
  int k = _width*i + j;
  assert(k < static_cast<int>(_S_cache.size()));
  if (_S_cache[k] < 0) {
    _S_cache[k] = _S(_h*j - _x0, _h*i - _y0);
  }
  return _S_cache[k];
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
