#ifndef __ABSTRACT_MARCHER_HPP__
#define __ABSTRACT_MARCHER_HPP__

struct abstract_marcher {
protected:
  abstract_marcher(double h, int S_cache_size);
  abstract_marcher(double h, double * S_cache);
  virtual ~abstract_marcher();
  double get_h() const { return _h; }
  double _h {1};
  double * _S_cache {nullptr};
private:
  bool _should_free_S_cache;
};

#endif // __ABSTRACT_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
