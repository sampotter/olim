#ifndef __COMMON_MACROS_HPP__
#define __COMMON_MACROS_HPP__

#define VAL(i) (nb[i]->get_value())
#define INF(x) std::numeric_limits<decltype(x)>::infinity()
#define ISINF(x) ((x) == INF(x))

#endif // __COMMON_MACROS_HPP__
