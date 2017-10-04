#ifndef __OLIM3D_MACROS_HPP__
#define __OLIM3D_MACROS_HPP__

#define RECT_LINE1(i) do {                        \
    T = min(T, this->line1(VAL(i), SPD(i), h));   \
  } while (0)

#define RECT_LINE2(i) do {                        \
    T = min(T, this->line2(VAL(i), SPD(i), h));   \
  } while (0)

#define RECT_LINE3(i) do {                        \
    T = min(T, this->line3(VAL(i), SPD(i), h));   \
  } while (0)

#define RECT_TRI11(i, j) do {                               \
    T = min(T, this->tri11(VAL(i), VAL(j), SPD(i, j), h));  \
  } while (0)

#define RECT_TRI12(i, j) do {                               \
    T = min(T, this->tri12(VAL(i), VAL(j), SPD(i, j), h));  \
  } while (0)

#define RECT_TRI13(i, j) do {                               \
    T = min(T, this->tri13(VAL(i), VAL(j), SPD(i, j), h));  \
  } while (0)

#define RECT_TRI22(i, j) do {                               \
    T = min(T, this->tri22(VAL(i), VAL(j), SPD(i, j), h));  \
  } while (0)

#define RECT_TRI23(i, j) do {                               \
    T = min(T, this->tri22(VAL(i), VAL(j), SPD(i, j), h));  \
  } while (0)

#define RECT_TETRA111(i, j, k) do {                                     \
    T = min(T, this->tetra111(VAL(i), VAL(j), VAL(k), SPD(i, j, k), h)); \
  } while (0)

#define RECT_TETRA122(i, j, k) do {                                     \
    T = min(T, this->tetra122(VAL(i), VAL(j), VAL(k), SPD(i, j, k), h)); \
  } while (0)

#define RECT_TETRA123(i, j, k) do {                                     \
    T = min(T, this->tetra123(VAL(i), VAL(j), VAL(k), SPD(i, j, k), h)); \
  } while (0)

#define RECT_TETRA222(i, j, k) do {                                     \
    T = min(T, this->tetra222(VAL(i), VAL(j), VAL(k), SPD(i, j, k), h)); \
  } while (0)

#endif // __OLIM3D_MACROS_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
