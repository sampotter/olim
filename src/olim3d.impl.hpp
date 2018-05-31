#ifndef __OLIM3D_IMPL_HPP__
#define __OLIM3D_IMPL_HPP__

#if COLLECT_STATS
#    include <cstdio>
#endif

#include <src/config.hpp>
#include "offsets.hpp"
#if COLLECT_STATS
#  include "update_rules.utils.hpp"
#endif

#define __di(l) di<3>[l]
#define __dj(l) dj<3>[l]
#define __dk(l) dk<3>[l]

/**
 * We only need to do triangle or line updates that lie in the xy-,
 * xz-, or yz-planes for every other octant using the preceding
 * enumeration of the octants (otherwise we double count these
 * updates). This function returns true if the passed octant requires
 * these updates.
 */
inline bool should_do_xyz_planar_updates(int octant) {
  static constexpr char parity[8] = {0, 1, 1, 0, 1, 0, 0, 1};
  return static_cast<bool>(parity[octant]);
}

namespace ind {
  // degree 1
  constexpr int N = 0;
  constexpr int E = 1;
  constexpr int U = 2;
  constexpr int S = 3;
  constexpr int W = 4;
  constexpr int D = 5;
  
  // degree 2
  constexpr int UN = 6;
  constexpr int UE = 7;
  constexpr int US = 8;
  constexpr int UW = 9;
  constexpr int NE = 10;
  constexpr int SE = 11;
  constexpr int SW = 12;
  constexpr int NW = 13;
  constexpr int DN = 14;
  constexpr int DE = 15;
  constexpr int DS = 16;
  constexpr int DW = 17;
  
  // degree 3
  constexpr int UNE = 18;
  constexpr int USE = 19;
  constexpr int USW = 20;
  constexpr int UNW = 21;
  constexpr int DNE = 22;
  constexpr int DSE = 23;
  constexpr int DSW = 24;
  constexpr int DNW = 25;
}

constexpr int oct2inds[8][7] = {
  {ind::D, ind::DS, ind::S, ind::SE, ind::E, ind::DE, ind::DSE},
  {ind::D, ind::DS, ind::S, ind::SW, ind::W, ind::DW, ind::DSW},
  {ind::D, ind::DN, ind::N, ind::NE, ind::E, ind::DE, ind::DNE},
  {ind::D, ind::DN, ind::N, ind::NW, ind::W, ind::DW, ind::DNW},
  {ind::U, ind::US, ind::S, ind::SE, ind::E, ind::UE, ind::USE},
  {ind::U, ind::US, ind::S, ind::SW, ind::W, ind::UW, ind::USW},
  {ind::U, ind::UN, ind::N, ind::NE, ind::E, ind::UE, ind::UNE},
  {ind::U, ind::UN, ind::N, ind::NW, ind::W, ind::UW, ind::UNW}
};

#define P001 1
#define P010 2
#define P011 3
#define P100 4
#define P101 5
#define P110 6
#define P111 7

// tri11 updates
#define __index02 0
#define __index20 0
#define __index24 1
#define __index42 1
#define __index04 2
#define __index40 2

// planar tri12 updates
#define __index01 3
#define __index10 3
#define __index05 4
#define __index50 4
#define __index21 5
#define __index12 5
#define __index23 6
#define __index32 6
#define __index34 7
#define __index43 7
#define __index45 8
#define __index54 8

// nonplanar tri12 updates
#define __index03 9
#define __index30 9
#define __index14 10
#define __index41 10
#define __index25 11
#define __index52 11

// tri13 updates
#define __index06 12
#define __index60 12
#define __index26 13
#define __index62 13
#define __index46 14
#define __index64 14

// tri22 updates
#define __index13 15
#define __index31 15
#define __index35 16
#define __index53 16
#define __index15 17
#define __index51 17

// tri23 updates
#define __index16 18
#define __index61 18
#define __index36 19
#define __index63 19
#define __index56 20
#define __index65 20

#define __is_planar_tri_update(i, j) ((__index##i##j) < 9)

#define __get_tri_skip_list(i, j)               \
  (__is_planar_tri_update(i, j) ? planar_tri_skip_list : tri_skip_list)

#define __get_tri_index(i, j)                                           \
  (__is_planar_tri_update(i, j) ? (__index##i##j) : (__index##i##j) - 9)

#define __get_offset(i, j) (__is_planar_tri_update(i, j) ? 9 : 12)

#define __skip_tri(i, j)                                    \
  __get_tri_skip_list(i, j)[__get_offset(i, j)*octant + __get_tri_index(i, j)]

#if COLLECT_STATS
#  define UPDATE_LINE_STATS(d) do {             \
    node_stats.inc_line_updates(d);             \
  } while (0)
#else
#  define UPDATE_LINE_STATS(d) do {} while (0)
#endif

#if COLLECT_STATS
#  define UPDATE_TRI_STATS(tmp, p0, p1) do {                            \
    tmp.hierarchical = l0 == parent || l1 == parent;                    \
    node_stats.inc_tri_updates(weight(p0), weight(p1),                  \
                               tmp.is_degenerate(),                     \
                               tmp.hierarchical);                       \
  } while (0)
#else
#  define UPDATE_TRI_STATS(tmp, p0, p1) do {} while (0)
#endif

#define TRI(i, j, p0, p1) do {                                  \
    if (!__skip_tri(i, j)) {                                    \
      int l0 = inds[i], l1 = inds[j];                           \
      if ((l0 == parent || l1 == parent) && nb[l0] && nb[l1]) { \
        auto tmp = this->tri(                                   \
          VAL(l0),                                              \
          VAL(l1),                                              \
          SPEED_ARGS(l0, l1),                                   \
          h,                                                    \
          ffvec<P##p0> {},                                      \
          ffvec<P##p1> {});                                     \
        T = min(T, tmp.value);                                  \
        UPDATE_TRI_STATS(tmp, P##p0, P##p1);                    \
      }                                                         \
    }                                                           \
  } while (0)

#if COLLECT_STATS
#  define UPDATE_TETRA_STATS(tmp, p0, p1, p2) do {                      \
    tmp.hierarchical = l0 == parent || l1 == parent || l2 == parent;    \
    node_stats.inc_tetra_updates(weight(p0), weight(p1), weight(p2),    \
                                 tmp.is_degenerate(),                   \
                                 tmp.hierarchical);                     \
  } while (0)
#else
#  define UPDATE_TETRA_STATS(tmp, p0, p1, p2) do {} while (0)
#endif

#define TETRA(i, j, k, p0, p1, p2) do {                                 \
    int l0 = inds[i], l1 = inds[j], l2 = inds[k];                       \
    if ((l0 == parent || l1 == parent || l2 == parent) &&               \
        nb[l0] && nb[l1] && nb[l2]) {                                   \
      auto tmp = this->tetra(                                           \
        VAL(l0),                                                        \
        VAL(l1),                                                        \
        VAL(l2),                                                        \
        SPEED_ARGS(l0, l1, l2),                                         \
        h,                                                              \
        ffvec<P ## p0> {},                                              \
        ffvec<P ## p1> {},                                              \
        ffvec<P ## p2> {});                                             \
      T = min(T, tmp.value);                                            \
      UPDATE_TETRA_STATS(tmp, P##p0, P##p1, P##p2);                     \
      __skip_tri(i, j) = 0x1;                                           \
      __skip_tri(j, k) = 0x1;                                           \
      __skip_tri(i, k) = 0x1;                                           \
    }                                                                   \
  } while (0)

#if COLLECT_STATS

template <class base_olim3d, class node, class line_updates, class tri_updates,
          class tetra_updates, int nneib>
void abstract_olim3d<
  base_olim3d, node, line_updates, tri_updates, tetra_updates,
  nneib>::dump_stats() const
{
  printf("depth = %d, width = %d, height = %d\n", this->get_depth(),
         this->get_width(), this->get_height());
  for (int k = 0; k < this->get_depth(); ++k) {
    for (int j = 0; j < this->get_width(); ++j) {
      for (int i = 0; i < this->get_height(); ++i) {
        auto const & stats = this->get_node_stats(i, j, k);
        printf("%d, %d, %d: visits = %d, line = %d, tri = %d/%d, tetra = %d/%d\n",
               i, j, k,
               stats.num_visits(),
               stats.num_line_updates(),
               stats.num_tri_updates() - stats.num_degenerate_tri_updates(),
               stats.num_tri_updates(),
               stats.num_tetra_updates() - stats.num_degenerate_tetra_updates(),
               stats.num_tetra_updates());
      }
    }
  }
}

#endif // COLLECT_STATS

template <class base_olim3d, class node, class line_updates, class tri_updates,
          class tetra_updates, int nneib>
void abstract_olim3d<
  base_olim3d, node, line_updates, tri_updates, tetra_updates,
  nneib>::get_valid_neighbors(int i, int j, int k, abstract_node ** nb)
{
  /**
   * TODO: conditionally only stage the neighbors that are necessary
   * based on which groups are being used
   */
  int a, b, c;
  for (int l = 0; l < nneib; ++l) {
    a = i + __di(l), b = j + __dj(l), c = k + __dk(l);
    if (this->in_bounds(a, b, c) && this->is_valid(a, b, c)) {
      nb[l] = &this->operator()(a, b, c);
    }
  }
}

template <
  class base_olim3d, class node, class line_updates, class tri_updates,
  class tetra_updates, int nneib>
void abstract_olim3d<
  base_olim3d, node, line_updates, tri_updates, tetra_updates,
  nneib>::stage_neighbors_impl(abstract_node * n)
{
  /**
   * TODO: conditionally only stage the neighbors that are necessary
   * based on which groups are being used
   */

  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();
  int k = static_cast<node *>(n)->get_k();
#if PRINT_UPDATES
  printf("olim3d::stage_neighbors_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  for (int l = 0; l < nneib; ++l) {
    this->stage(i + __di(l), j + __dj(l), k + __dk(l));
  }

  int a, b, c, l;
  for (l = 0; l < std::min(6, nneib); ++l) {
    a = i + __di(l), b = j + __dj(l), c = k + __dk(l);
    if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      this->update(a, b, c, (l + 3) % 6);
    }
  }
  for (l = 6; l < std::min(18, nneib); ++l) {
    a = i + __di(l), b = j + __dj(l), c = k + __dk(l);
    if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      this->update(a, b, c, 22 - 2*(l/2) + (l % 2));
    }
  }
  for (l = 18; l < std::min(26, nneib); ++l) {
    a = i + __di(l), b = j + __dj(l), c = k + __dk(l);
    if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      this->update(a, b, c, 42 - 2*(l/2) + (l % 2));
    }
  }
}

template <
  class base_olim3d, class node, class line_updates, class tri_updates,
  class tetra_updates, int nneib>
void abstract_olim3d<
  base_olim3d, node, line_updates, tri_updates, tetra_updates,
  nneib>::update_impl(int i, int j, int k, int parent, double & T)
{
  static_cast<base_olim3d *>(this)->update_crtp(i, j, k, parent, T);
}

template <class node, class line_updates, class tri_updates,
          class tetra_updates, class groups>
void olim3d_bv<
  node, line_updates, tri_updates, tetra_updates,
  groups>::update_crtp(int i, int j, int k, int parent, double & T)
{
  // TODO: not currently using this. An easy way to use it would be to
  // map each parent index to a list of octants to iterate over.
  (void) parent;

  using std::min;
#if PRINT_UPDATES
  printf("olim3d::update_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  abstract_node * nb[groups::nneib];
  memset(nb, 0x0, groups::nneib*sizeof(abstract_node *));
  this->get_valid_neighbors(i, j, k, nb);

  // TODO: once we know the number of neighbors that are actually
  // valid, we can alloca and reduce the size of the arrays that
  // follow (just create an array of indices)
  double h = this->get_h(), s = this->get_speed(i, j, k), s_[groups::nneib];
  for (int l = 0; l < groups::nneib; ++l) {
    if (nb[l]) {
      s_[l] = this->get_speed(i + __di(l), j + __dj(l), k + __dk(l));
    }
  }

#if COLLECT_STATS
  auto & node_stats = this->get_node_stats(i, j, k);
  node_stats.inc_num_visits();
#endif

  // Do line update corresponding to parent node.
  if (parent < 6) {
    T = min(T, this->template line<1>(VAL(parent), SPEED_ARGS(parent), h));
    UPDATE_LINE_STATS(1);
  } else if (parent < 18) {
    T = min(T, this->template line<2>(VAL(parent), SPEED_ARGS(parent), h));
    UPDATE_LINE_STATS(2);
  } else {
    T = min(T, this->template line<3>(VAL(parent), SPEED_ARGS(parent), h));
    UPDATE_LINE_STATS(3);
  }

  int const * inds;

  /**
   * These arrays are used to keep track of which triangle and line
   * updates should be skipped.
   */
  char planar_tri_skip_list[8*9]; // this is inefficient--only need [4][9]
  char tri_skip_list[8*12];
  memset((void *) planar_tri_skip_list, 0x0, sizeof(char)*8*9);
  memset((void *) tri_skip_list, 0x0, sizeof(char)*8*12);

  /**
   * Tetrahedron updates:
   */
  for (int octant = 0; octant < 8; ++octant) {
    inds = oct2inds[octant];
    if (groups::group_I) {
      TETRA(1, 2, 3, 011, 010, 110);
      TETRA(3, 4, 5, 110, 100, 101);
      TETRA(5, 0, 1, 101, 001, 011);
    }
    if (groups::group_II) {
      TETRA(0, 1, 3, 001, 011, 110);
      TETRA(1, 2, 4, 011, 010, 100);
      TETRA(2, 3, 5, 010, 110, 101);
      TETRA(3, 4, 0, 110, 100, 001);
      TETRA(4, 5, 1, 100, 101, 011);
      TETRA(5, 0, 2, 101, 001, 010);
    }
    if (groups::group_III) {
      TETRA(0, 1, 4, 001, 011, 100);
      TETRA(1, 2, 5, 011, 010, 101);
      TETRA(2, 3, 0, 010, 110, 001);
      TETRA(3, 4, 1, 110, 100, 011);
      TETRA(4, 5, 2, 100, 101, 010);
      TETRA(5, 0, 3, 101, 001, 110);
    }
    if (groups::group_IV_a) {
      TETRA(0, 2, 4, 001, 010, 100);
    }
    if (groups::group_IV_b) {
      TETRA(1, 3, 5, 011, 110, 101);
    }
    if (groups::group_V) {
      TETRA(0, 1, 6, 001, 011, 111);
      TETRA(1, 2, 6, 011, 010, 111);
      TETRA(2, 3, 6, 010, 110, 111);
      TETRA(3, 4, 6, 110, 100, 111);
      TETRA(4, 5, 6, 100, 101, 111);
      TETRA(5, 0, 6, 101, 001, 111);
    }
    if (groups::group_VI_a) {
      TETRA(0, 2, 6, 001, 010, 111);
      TETRA(2, 4, 6, 010, 100, 111);
      TETRA(4, 0, 6, 100, 001, 111);
    }
    if (groups::group_VI_b) {
      TETRA(1, 3, 6, 011, 110, 111);
      TETRA(3, 5, 6, 110, 101, 111);
      TETRA(5, 1, 6, 101, 011, 111);
    }
  }

  /**
   * Triangle updates:
   */
  for (int octant = 0; octant < 8; ++octant) {
    inds = oct2inds[octant];
    if (should_do_xyz_planar_updates(octant)) {
      if (groups::do_tri11_updates) {
        TRI(0, 2, 001, 010);
        TRI(2, 4, 010, 100);
        TRI(4, 0, 100, 001);
      }
      if (groups::do_tri12_updates) {
        TRI(0, 1, 001, 011);
        TRI(2, 1, 010, 011);
        TRI(2, 3, 010, 110);
        TRI(4, 3, 100, 110);
        TRI(4, 5, 100, 101);
        TRI(0, 5, 001, 101);
      }
    }
    if (groups::do_tri13_updates) {
      TRI(0, 6, 001, 111);
      TRI(2, 6, 010, 111);
      TRI(4, 6, 100, 111);
    }
    if (groups::do_tri22_updates) {
      TRI(1, 3, 011, 110);
      TRI(3, 5, 110, 101);
      TRI(5, 1, 101, 011);
    }
    if (groups::do_tri23_updates) {
      TRI(1, 6, 011, 111);
      TRI(3, 6, 110, 111);
      TRI(5, 6, 101, 111);
    }
  }

#if PRINT_UPDATES
  printf("olim3d::update_impl: T <- %g\n", T);
#endif
}

template <class node, class line_updates, class tri_updates,
          class tetra_updates, int lp_norm, int d1, int d2>
void olim3d_hu<
  node, line_updates, tri_updates, tetra_updates,
  lp_norm, d1, d2>::update_crtp(int i, int j, int k, int parent, double & T)
{
  using std::min;
#if PRINT_UPDATES
  printf("olim3d_hu::update_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  // TODO: an idea for an optimization: count the number of valid
  // neighbors, and create an array of indices of that size on the
  // stack using `alloca' storing the indices of the valid
  // neighbors. Make s[] the same size. Then, we only have to deal
  // with the valid neighbors directly and don't need to check "if
  // (nb[l]) { ... }" over and over again.

  abstract_node * nb[26];
  memset(nb, 0x0, 26*sizeof(abstract_node *));
  this->get_valid_neighbors(i, j, k, nb);

  double h = this->get_h(), s = this->get_speed(i, j, k), s_[26];
  for (int l = 0; l < 26; ++l) {
    if (nb[l]) {
      s_[l] = this->get_speed(i + __di(l), j + __dj(l), k + __dk(l));
    }
  }

#if COLLECT_STATS
  auto & node_stats = this->get_node_stats(i, j, k);
  node_stats.inc_num_visits();
#endif

  /**
   * Tnew: temporary variable used for updating T, the update value
   * T0, T1, T2: separate minimum values for degree 0/1/2 updates
   * l0, l1: indices of p0 and p1
   * p0, p1, p2: update node vectors
   */
  double Tnew, T0 = INF(double), T1 = INF(double), T2 = INF(double);
  int l0 = -1, l1 = -1;
  double p0[3], p1[3], p2[3];

  // Depending on the compilation flag HU_USE_PARENT_NODE, either let
  // l0 be the index of the updating parent, or let l0 be the index of
  // the neighboring valid node with the minimum degree 0 update
  // value.
#if HU_USE_PARENT_NODE
  l0 = parent;
  p0[0] = __di(l0);
  p0[1] = __dj(l0);
  p0[2] = __dk(l0);
  T0 = this->template line<3>(p0, VAL(l0), SPEED_ARGS(l0), h);
#  if COLLECT_STATS
  node_stats.inc_line_updates(p0, 3);
#  endif
#else // HU_USE_PARENT_NODE
  for (int l = 0; l < 26; ++l) {
    if (nb[l]) {
      p0[0] = __di(l);
      p0[1] = __dj(l);
      p0[2] = __dk(l);
      Tnew = this->template line<3>(p0, VAL(l), SPEED_ARGS(l), h);
#  if COLLECT_STATS
      node_stats.inc_line_updates(p0, 3);
#  endif
      if (Tnew < T0) {
        T0 = Tnew;
        l0 = l;
      }
    }
  }
  assert(l0 != -1);
  p0[0] = __di(l0);
  p0[1] = __dj(l0);
  p0[2] = __dk(l0);
#endif // HU_USE_PARENT_NODE

  // Create a cache for the minimizing lambdas to use for skipping
  // tetrahedron updates. Initialize it to -1 so that we can tell
  // which triangle updates have been computed.
  double arglam[26];
  std::fill(arglam, arglam + sizeof(arglam)/sizeof(double), -1);

  // Find the minimal triangle update containing l0.
  for (int l = 0; l < 26; ++l) {
    if (l != l0 && nb[l]) {
      p1[0] = __di(l);
      p1[1] = __dj(l);
      p1[2] = __dk(l);

      // Check to see how far apart p0 and p1, and continue to the
      // next candidate for p1 if they're too far apart.
      if (lp_norm == L1) {
        if (dist1<3>(p0, p1) > d1) continue;
      } else if (lp_norm == L2) {
        if (dist2sq<3>(p0, p1) > d1) continue;
      } else {
        if (distmax<3>(p0, p1) > d1) continue;
      }

      // Do the triangle update.
      auto const tmp = this->template tri<3>(
        p0, p1, VAL(l0), VAL(l), SPEED_ARGS(l0, l), h);
#if COLLECT_STATS
      node_stats.inc_tri_updates(p0, p1, 3, tmp.is_degenerate(), true);
#endif
      Tnew = tmp.value;
      if (Tnew < T1) {
        T1 = Tnew;
        l1 = l;
      }
      arglam[l] = tmp.lambda[0];
    }
  }
  // There may only be a single valid neighbor, in which case we
  // should jump to the end of this function and skip all of the
  // tetrahedron updates.
  if (l1 == -1) {
    assert(std::isinf(T1));
    goto coda;
  }
  p1[0] = __di(l1);
  p1[1] = __dj(l1);
  p1[2] = __dk(l1);

  {
    // Use the scalar triple product (dot(p x q, r)) to check if three
    // points are coplanar. We precompute the cross product of p0 and
    // p1 here so that we only have to compute a dot product up ahead.
    double const p0_cross_p1[3] = {
      p0[1]*p1[2] - p0[2]*p1[1],
      p0[2]*p1[0] - p0[0]*p1[2],
      p0[0]*p1[1] - p0[1]*p1[0]
    };

    // Do the tetrahedron updates such that p2 is sufficiently near p0
    // and p1.
    for (int l2 = 0; l2 < 26; ++l2) {
      if (l2 != l1 && l2 != l0 && nb[l2]) {
        p2[0] = __di(l2);
        p2[1] = __dj(l2);
        p2[2] = __dk(l2);

        // Check if p0, p1, and p2 are coplanar: if they are, we can
        // skip the tetrahedron update.
        if (fabs(p0_cross_p1[0]*p2[0] +
                 p0_cross_p1[1]*p2[1] +
                 p0_cross_p1[2]*p2[2]) < 1e2*EPS(double)) {
          continue;
        }

        // Check to see how far p2 is from p0 and p1---if it's too
        // far, proceed to the next candidate for p2.
        if (lp_norm == L1) {
          if (dist1<3>(p0, p2) > d2 || dist1<3>(p1, p2) > d2) continue;
        } else if (lp_norm == L2) {
          if (dist2sq<3>(p0, p2) > d2 || dist2sq<3>(p1, p2) > d2) continue;
        } else {
          if (distmax<3>(p0, p2) > d2 || distmax<3>(p1, p2) > d2) continue;
        }

        // Compute Lagrange multipliers to see if we can skip the
        // tetrahedron update.

        // TODO: this is pretty inefficient! we will just end up
        // redoing the setup phase for this cost_func once we get
        // inside tetra in the next section
        typename tetra_updates::template cost_func<3, 2> func(h, this->theta());
        {
          double U[3] = {VAL(l0), VAL(l1), VAL(l2)};
          double S[3] = {s_[l0], s_[l1], s_[l2]};
          double P[3][3] = {
            {p0[0], p0[1], p0[2]},
            {p1[0], p1[1], p1[2]},
            {p2[0], p2[1], p2[2]}
          };
          func.set_args(U, s, S, P);
        }

        double mu[2], lam[2];
        int k;

        // We get the first two checks for free using arglam

        assert(arglam[l1] != -1);
        lam[0] = arglam[l1];
        lam[1] = 0;
        func.lag_mult(lam, mu, &k);
        if (mu[0] < 0 || (k == 2 && mu[1] < 0)) {
          continue;
        }

        assert(arglam[l2] != -1);
        lam[0] = 0;
        lam[1] = arglam[l2];
        func.lag_mult(lam, mu, &k);
        if (mu[0] < 0 || (k == 2 && mu[1] < 0)) {
          continue;
        }

        // TODO: We're not doing the third triangle update right
        // now. This seems to work okay for now, but we can't do the
        // "exact solve" using the QR decomposition if we don't do the
        // third update. It's unclear how efficient this will be...

        // Finally, do the tetrahedron update.
        auto const tmp = this->template tetra(
          p0, p1, p2, VAL(l0), VAL(l1), VAL(l2), SPEED_ARGS(l0, l1, l2), h);
#if COLLECT_STATS
        node_stats.inc_tetra_updates(p0, p1, p2, 3, tmp.is_degenerate(), true);
#endif
        Tnew = tmp.value;
        if (Tnew < T2) {
          T2 = Tnew;
        }
      }
    }
  }

  // Set T to be the minimum of T0, T1, and T2.
coda:
  T = min(T0, min(T1, T2));

#if PRINT_UPDATES
  printf("olim3d_hu::update_impl: T <- %g\n", T);
#endif
}

#undef LINE
#undef TRI
#undef TETRA

#undef __is_planar_tri_update
#undef __get_tri_skip_list
#undef __get_tri_index
#undef __get_offset
#undef __skip_tri

#undef P001
#undef P010
#undef P011
#undef P100
#undef P101
#undef P110
#undef P111

#undef __di
#undef __dj
#undef __dk

#if COLLECT_STATS

template <class node, class line_updates, class tri_updates,
          class tetra_updates, class groups, int nneib>
olim3d_node_stats &
abstract_olim3d<
  node, line_updates, tri_updates, tetra_updates, groups,
  nneib>::get_node_stats(int i, int j, int k)
{
  return _node_stats[this->get_height()*(this->get_width()*k + j) + i];
}

template <class node, class line_updates, class tri_updates,
          class tetra_updates, class groups, int nneib>
olim3d_node_stats const &
abstract_olim3d<
  node, line_updates, tri_updates, tetra_updates, groups,
  nneib>::get_node_stats(int i, int j, int k) const
{
  return _node_stats[this->get_height()*(this->get_width()*k + j) + i];
}

#endif // COLLECT_STATS

#endif // __OLIM3D_IMPL_HPP__
