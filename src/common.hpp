#ifndef __COMMON_HPP__
#define __COMMON_HPP__

namespace eikonal {
  template <bool value>
  using bool_t = std::integral_constant<bool, value>;
}

#endif // __COMMON_HPP__
