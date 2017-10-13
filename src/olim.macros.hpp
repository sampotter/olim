#ifndef __OLIM_MACROS_HPP__
#define __OLIM_MACROS_HPP__

#include <src/config.hpp>

#define LINE1(k) do {                                                  \
    T = min(T, this->line1(VAL(k), s, speed(i + di[k], j + dj[k]), h)); \
  } while (0)

#define LINE2(k) do {                                                   \
    T = min(T, this->line2(VAL(k), s, speed(i + di[k], j + dj[k]), h)); \
  } while (0)                                                           \

#ifdef OLIM8_ADJ_UDPATES
#define TRI12(k, l) do {                                    \
    T = min(T, this->tri12(VAL(k), VAL(l), s,               \
                           speed(i + di[k], j + dj[k]),     \
                           speed(i + di[l], j + dj[l]),     \
                           h));                             \
  } while (0)
#endif

#define TRI22(k, l) do {                                    \
    T = min(T, this->tri22(VAL(k), VAL(l), s,               \
                           speed(i + di[k], j + dj[k]),     \
                           speed(i + di[l], j + dj[l]),     \
                           h));                             \
  } while (0)

#endif // __OLIM_MACROS_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
