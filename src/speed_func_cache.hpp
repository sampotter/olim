#ifndef __SPEED_FUNC_CACHE_HPP__
#define __SPEED_FUNC_CACHE_HPP__

struct speed_func_cache {
  speed_func_cache(int size, double * cache_ptr = nullptr);
  ~speed_func_cache();

  inline bool is_cached(int i) const { return _cache[i] >= 0; }
  inline double cached_value(int i) const { return _cache[i]; }
  inline void cache_value(int i, double x) { _cache[i] = x; }
private:
  bool _should_free;
  double * _cache;
};

#endif // __SPEED_FUNC_CACHE_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
