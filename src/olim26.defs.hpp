#ifndef __OLIM26_DEFS_HPP_HPP__
#define __OLIM26_DEFS_HPP_HPP__

namespace olim26_defs {
  enum dir {
    N, E, U, S, W, D, // degree 1
 // 0  1  2  3  4  5
    UN, UE, US, UW, NE, SE, SW, NW, DN, DE, DS, DW, // degree 2
 // 6   7   8   9   10  11  12  13  14  15  16  17
    UNE, USE, USW, UNW, DNE, DSE, DSW, DNW, // degree 3
 // 18   19   20   21   22   23   24   25
  };
}

#endif // __OLIM26_DEFS_HPP_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
