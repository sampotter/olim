#ifndef __COMMON_MACROS_HPP__
#define __COMMON_MACROS_HPP__

#define EPS(T) std::numeric_limits<T>::epsilon()
#define INF(T) std::numeric_limits<T>::infinity()
#define VAL(i) (nb[i]->get_value())

#define GET_MACRO_NAME_3(_1, _2, _3, NAME, ...) NAME

#endif // __COMMON_MACROS_HPP__
