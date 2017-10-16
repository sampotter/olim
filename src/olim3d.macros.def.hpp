#ifdef OLIM_MACROS_DEFINED
#    error "OLIM macros have already been defined"
#endif

#define OLIM_MACROS_DEFINED 1

#define LINE1(i) do {                                   \
    T = min(T, this->line1(VAL(i), SPEED_ARGS(i), h));  \
  } while (0)

#define LINE2(i) do {                                   \
    T = min(T, this->line2(VAL(i), SPEED_ARGS(i), h));  \
  } while (0)

#define LINE3(i) do {                                   \
    T = min(T, this->line3(VAL(i), SPEED_ARGS(i), h));  \
  } while (0)

#define TRI11(i, j) do {                                    \
    T = min(                                                \
      T, this->tri11(                                       \
        VAL(i),                                             \
        VAL(j),                                             \
        SPEED_ARGS(i, j),                                   \
        h));                                                \
  } while (0)

#define TRI12(i, j) do {                                    \
    T = min(                                                \
      T, this->tri12(                                       \
        VAL(i),                                             \
        VAL(j),                                             \
        SPEED_ARGS(i, j),                                   \
        h));                                                \
  } while (0)

#define TRI13(i, j) do {                                    \
    T = min(                                                \
      T, this->tri13(                                       \
        VAL(i),                                             \
        VAL(j),                                             \
        SPEED_ARGS(i, j),                                   \
        h));                                                \
  } while (0)

#define TRI22(i, j) do {                                    \
    T = min(                                                \
      T, this->tri22(                                       \
        VAL(i),                                             \
        VAL(j),                                             \
        SPEED_ARGS(i, j),                                   \
        h));                                                \
  } while (0)

#define TRI23(i, j) do {                                    \
    T = min(                                                \
      T, this->tri23(                                       \
        VAL(i),                                             \
        VAL(j),                                             \
        SPEED_ARGS(i, j),                                   \
        h));                                                \
  } while (0)

#define TETRA111(i, j, k) do {                                          \
    T = min(                                                            \
      T,                                                                \
      this->tetra111(                                                   \
        VAL(i),                                                         \
        VAL(j),                                                         \
        VAL(k),                                                         \
        SPEED_ARGS(i, j, k),                                            \
        h));                                                            \
  } while (0)

#define TETRA122(i, j, k) do {                                          \
    T = min(                                                            \
      T,                                                                \
      this->tetra122(                                                   \
        VAL(i),                                                         \
        VAL(j),                                                         \
        VAL(k),                                                         \
        SPEED_ARGS(i, j, k),                                            \
        h));                                                            \
  } while (0)

#define TETRA123(i, j, k) do {                                          \
    T = min(                                                            \
      T,                                                                \
      this->tetra123(                                                   \
        VAL(i),                                                         \
        VAL(j),                                                         \
        VAL(k),                                                         \
        SPEED_ARGS(i, j, k),                                            \
        h));                                                            \
  } while (0)

#define TETRA222(i, j, k) do {                                          \
    T = min(                                                            \
      T,                                                                \
      this->tetra222(                                                   \
        VAL(i),                                                         \
        VAL(j),                                                         \
        VAL(k),                                                         \
        SPEED_ARGS(i, j, k),                                            \
        h));                                                            \
  } while (0)

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
