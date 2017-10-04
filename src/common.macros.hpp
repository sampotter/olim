#ifndef __COMMON_MACROS_HPP__
#define __COMMON_MACROS_HPP__

#define EPS(T) std::numeric_limits<T>::epsilon()
#define INF(T) std::numeric_limits<T>::infinity()
#define VAL(i) (nb[i]->get_value())

#define GET_MACRO_NAME_3(_1, _2, _3, NAME, ...) NAME
#define SPD(...) GET_MACRO_NAME_3(__VA_ARGS__, SPD3, SPD2, SPD1)(__VA_ARGS__)
#define SPD1(i) (this->estimate_speed(s, s_[i]))
#define SPD2(i, j) (this->estimate_speed(s, s_[i], s_[j]))
#define SPD3(i, j, k) (this->estimate_speed(s, s_[i], s_[j], s_[k]))

#endif // __COMMON_MACROS_HPP__
