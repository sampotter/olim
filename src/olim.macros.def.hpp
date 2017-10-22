#ifdef OLIM_MACROS_DEFINED
#    error "OLIM macros have already been defined"
#endif

#define OLIM_MACROS_DEFINED 1

#include <src/config.hpp>

#define LINE1(k) do {                           \
    T = min(                                    \
      T,                                        \
      this->template line<1>(                   \
        VAL(k),                                 \
        s,                                      \
        speed(i + di[k], j + dj[k]),            \
        h));                                    \
  } while (0)

#define LINE2(k) do {                           \
    T = min(                                    \
      T,                                        \
      this->template line<2>(                   \
        VAL(k),                                 \
        s,                                      \
        speed(i + di[k], j + dj[k]),            \
        h));                                    \
  } while (0)                                   \

#ifdef OLIM8_ADJ_UDPATES
#define TRI11(k, l) do {                                     \
    T = min(T, this->adj2pt(VAL(k), VAL(l), s,               \
                            speed(i + di[k], j + dj[k]),     \
                            speed(i + di[l], j + dj[l]),     \
                            h));                             \
  } while (0)
#endif

#define TRI12(k, l) do {                                      \
    T = min(T, this->tri12(VAL(k), VAL(l), s,                 \
                           speed(i + di[k], j + dj[k]),       \
                           speed(i + di[l], j + dj[l]),       \
                           h));                               \
  } while (0)

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
