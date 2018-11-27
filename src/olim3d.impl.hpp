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

// // tri11 updates
// #define __index02 0
// #define __index20 0
// #define __index24 1
// #define __index42 1
// #define __index04 2
// #define __index40 2

// // tri12 updates
// #define __index01 3
// #define __index10 3
// #define __index05 4
// #define __index50 4
// #define __index21 5
// #define __index12 5
// #define __index23 6
// #define __index32 6
// #define __index34 7
// #define __index43 7
// #define __index45 8
// #define __index54 8
// #define __index03 9
// #define __index30 9
// #define __index14 10
// #define __index41 10
// #define __index25 11
// #define __index52 11

// // tri13 updates
// #define __index06 12
// #define __index60 12
// #define __index26 13
// #define __index62 13
// #define __index46 14
// #define __index64 14

// // tri22 updates
// #define __index13 15
// #define __index31 15
// #define __index35 16
// #define __index53 16
// #define __index15 17
// #define __index51 17

// // tri23 updates
// #define __index16 18
// #define __index61 18
// #define __index36 19
// #define __index63 19
// #define __index56 20
// #define __index65 20

// #define __skip_tri(i, j) tri_skip_list[21*octant + __index##i##j]

// #if COLLECT_STATS
// #  define UPDATE_TRI_STATS(tmp, p0, p1) do {                            \
//     tmp.hierarchical = l0 == parent || l1 == parent;                    \
//     node_stats.inc_tri_updates(weight(p0), weight(p1),                  \
//                                tmp.is_degenerate(),                     \
//                                tmp.hierarchical);                       \
//   } while (0)
// #else
// #  define UPDATE_TRI_STATS(tmp, p0, p1) do {} while (0)
// #endif

// #if TRACK_PARENTS
// #  define __tri_update_value() do {             \
//     if (tmp.value < T) {                        \
//       T = tmp.value;                            \
//       n->set_parents({{                         \
//         &this->operator()(                      \
//           n->get_i() + __di(l0),                \
//           n->get_j() + __dj(l0),                \
//           n->get_k() + __dk(l0)),               \
//         &this->operator()(                      \
//           n->get_i() + __di(l1),                \
//           n->get_j() + __dj(l1),                \
//           n->get_k() + __dk(l1)),               \
//         nullptr                                 \
//       }});                                      \
//     }                                           \
//   } while (0)
// #else
// #  define __tri_update_value() do {             \
//     T = min(T, tmp.value);                      \
//   } while (0)
// #endif

// #define TRI(i, j, p0, p1) do {                                  \
//     if (!__skip_tri(i, j)) {                                    \
//       int l0 = inds[i], l1 = inds[j];                           \
//       if ((l0 == parent || l1 == parent) && nb[l0] && nb[l1]) { \
//         auto tmp = this->tri(                                   \
//           VAL(l0),                                              \
//           VAL(l1),                                              \
//           SPEED_ARGS(l0, l1),                                   \
//           h,                                                    \
//           ffvec<P##p0> {},                                      \
//           ffvec<P##p1> {});                                     \
//         __tri_update_value();                                   \
//         UPDATE_TRI_STATS(tmp, P##p0, P##p1);                    \
//       }                                                         \
//       __skip_tri(i, j) = 0x1;                                   \
//     }                                                           \
//   } while (0)

// #define TRI_FAC(i, j) do {                                      \
//     if (!__skip_tri(i, j)) {                                    \
//       int l0 = inds[i], l1 = inds[j];                           \
//       if ((l0 == parent || l1 == parent) && nb[l0] && nb[l1]) { \
//         p0[0] = __di(l0);                                       \
//         p0[1] = __dj(l0);                                       \
//         p0[2] = __dk(l0);                                       \
//         p1[0] = __di(l1);                                       \
//         p1[1] = __dj(l1);                                       \
//         p1[2] = __dk(l1);                                       \
//         auto tmp = this->template tri<3>(                       \
//           VAL(l0),                                              \
//           VAL(l1),                                              \
//           SPEED_ARGS(l0, l1),                                   \
//           h,                                                    \
//           p0,                                                   \
//           p1,                                                   \
//           p_fac,                                                \
//           s_fac);                                               \
//         __tri_update_value();                                   \
//         /* TODO: collect stats */                               \
//       }                                                         \
//       __skip_tri(i, j) = 0x1;                                   \
//     }                                                           \
//   } while (0)

// #if COLLECT_STATS
// #  define UPDATE_TETRA_STATS(tmp, p0, p1, p2) do {                      \
//     tmp.hierarchical = l0 == parent || l1 == parent || l2 == parent;    \
//     node_stats.inc_tetra_updates(weight(p0), weight(p1), weight(p2),    \
//                                  tmp.is_degenerate(),                   \
//                                  tmp.hierarchical);                     \
//   } while (0)
// #else
// #  define UPDATE_TETRA_STATS(tmp, p0, p1, p2) do {} while (0)
// #endif

// #if TRACK_PARENTS
// #  define __tetra_update_value() do {                                   \
//     if (tmp.value < T) {                                                \
//       T = tmp.value;                                                    \
//       n->set_parents({{                                                 \
//         &this->operator()(                                              \
//           n->get_i() + __di(l0),                                        \
//           n->get_j() + __dj(l0),                                        \
//           n->get_k() + __dk(l0)),                                       \
//         &this->operator()(                                              \
//           n->get_i() + __di(l1),                                        \
//           n->get_j() + __dj(l1),                                        \
//           n->get_k() + __dk(l1)),                                       \
//         &this->operator()(                                              \
//           n->get_i() + __di(l2),                                        \
//           n->get_j() + __dj(l2),                                        \
//           n->get_k() + __dk(l2))                                        \
//       }});                                                              \
//     }                                                                   \
//   } while (0)
// #else
// #  define __tetra_update_value() do {           \
//     T = min(T, tmp.value);                      \
//   } while (0)
// #endif

// #define TETRA(i, j, k, p0, p1, p2) do {                                 \
//     int l0 = inds[i], l1 = inds[j], l2 = inds[k];                       \
//     if ((l0 == parent || l1 == parent || l2 == parent) &&               \
//         nb[l0] && nb[l1] && nb[l2]) {                                   \
//       auto tmp = this->tetra(                                           \
//         VAL(l0),                                                        \
//         VAL(l1),                                                        \
//         VAL(l2),                                                        \
//         SPEED_ARGS(l0, l1, l2),                                         \
//         h,                                                              \
//         ffvec<P ## p0> {},                                              \
//         ffvec<P ## p1> {},                                              \
//         ffvec<P ## p2> {});                                             \
//       __tetra_update_value();                                           \
//       UPDATE_TETRA_STATS(tmp, P##p0, P##p1, P##p2);                     \
//       __skip_tri(i, j) = 0x1;                                           \
//       __skip_tri(j, k) = 0x1;                                           \
//       __skip_tri(i, k) = 0x1;                                           \
//     }                                                                   \
//   } while (0)

// #define TETRA_FAC(i, j, k) do {                                         \
//     int l0 = inds[i], l1 = inds[j], l2 = inds[k];                       \
//     if ((l0 == parent || l1 == parent || l2 == parent) &&               \
//         nb[l0] && nb[l1] && nb[l2]) {                                   \
//       p0[0] = __di(l0);                                                 \
//       p0[1] = __dj(l0);                                                 \
//       p0[2] = __dk(l0);                                                 \
//       p1[0] = __di(l1);                                                 \
//       p1[1] = __dj(l1);                                                 \
//       p1[2] = __dk(l1);                                                 \
//       p2[0] = __di(l2);                                                 \
//       p2[1] = __dj(l2);                                                 \
//       p2[2] = __dk(l2);                                                 \
//       auto tmp = this->tetra(                                           \
//         VAL(l0),                                                        \
//         VAL(l1),                                                        \
//         VAL(l2),                                                        \
//         SPEED_ARGS(l0, l1, l2),                                         \
//         h,                                                              \
//         p0,                                                             \
//         p1,                                                             \
//         p2,                                                             \
//         p_fac,                                                          \
//         s_fac);                                                         \
//       __tetra_update_value();                                           \
//       /* TODO: collect stats */                                         \
//       __skip_tri(i, j) = 0x1;                                           \
//       __skip_tri(j, k) = 0x1;                                           \
//       __skip_tri(i, k) = 0x1;                                           \
//     }                                                                   \
//   } while (0)

#if COLLECT_STATS

template <class base, class node, int num_neighbors>
void abstract_olim3d<base, node, num_neighbors>::dump_stats() const
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

template <class base, class node, int num_neighbors>
void abstract_olim3d<base, node, num_neighbors>::update_impl(
  node * n, node ** nb, int parent, double & T)
{
  int i = n->get_i(), j = n->get_j(), k = n->get_k();
  for (int l = 0; l < num_neighbors; ++l) {
    if (nb[l]) {
      this->s[l] = this->get_speed(i + __di(l), j + __dj(l), k + __dk(l));
    }
  }

  static_cast<base *>(this)->n = n;
  static_cast<base *>(this)->nb = nb;
  static_cast<base *>(this)->parent = parent;
  static_cast<base *>(this)->update_crtp(T);
}

template <cost_func F, class node, class groups>
void olim3d_bv<F, node, groups>::update_crtp(double & T)
{
  using std::min;

  assert(parent <= groups::num_neighbors);

  // int i = n->get_i(), j = n->get_j(), k = n->get_k();
#if PRINT_UPDATES
  printf("olim3d::update_impl(i = %d, j = %d, k = %d)\n",
         n->get_i(), n->get_j(), n->get_k());
#endif

#if COLLECT_STATS
  node_stats = this->get_node_stats(i, j, k);
  node_stats.inc_num_visits();
#endif

  // Do line update corresponding to parent node.
  {
    double Tnew = INF(double);
    if (parent < 6) {
      line<1>(parent, Tnew);
    } else if (groups::num_neighbors > 6 && parent < 18) {
      line<2>(parent, Tnew);
    } else if (groups::num_neighbors > 18) {
      line<3>(parent, Tnew);
    }
    assert(!std::isinf(Tnew));
#if TRACK_PARENTS
    if (Tnew < T) {
      T = Tnew;
      n->set_parents({{
        &this->operator()(i + __di(parent), j + __dj(parent), k + __dk(parent)),
        nullptr,
        nullptr}});
    }
#else
    T = min(T, Tnew);
#endif
  }

  /**
   * This array is used to keep track of which triangle updates should
   * be skipped.
   *
   * TODO: we could use a bitvector here instead.
   */
  // char tri_skip_list[8*21];
  // memset((void *) tri_skip_list, 0x0, sizeof(char)*8*21);
  reset_tri_skip_list();

  if (n->has_fac_parent()) {
    // auto n_fac = static_cast<node *>(n->get_fac_parent());
    // int i_fac = n_fac->get_i(), j_fac = n_fac->get_j(), k_fac = n_fac->get_k();
    // double s_fac = this->get_speed(i_fac, j_fac, k_fac);
    // // double p0[3], p1[3], p2[3]; // TODO: these should be template parameters!
    // double p_fac[3] = {
    //   (double) (i_fac - i),
    //   (double) (j_fac - j),
    //   (double) (k_fac - k)
    // };

    /**
     * Tetrahedron updates:
     */
    // TODO: decide which octants to do based on the parent index
    // (only need to look at 1, 2, or 4 octants, never more than that)
    for (octant = 0; octant < 8; ++octant) {
      inds = oct2inds[octant];
      if (groups::group_I) {
        tetra_fac<1, 2, 3>(T);
        tetra_fac<3, 4, 5>(T);
        tetra_fac<5, 0, 1>(T);
      }
      if (groups::group_II) {
        tetra_fac<0, 1, 3>(T);
        tetra_fac<1, 2, 4>(T);
        tetra_fac<2, 3, 5>(T);
        tetra_fac<3, 4, 0>(T);
        tetra_fac<4, 5, 1>(T);
        tetra_fac<5, 0, 2>(T);
      }
      if (groups::group_III) {
        tetra_fac<0, 1, 4>(T);
        tetra_fac<1, 2, 5>(T);
        tetra_fac<2, 3, 0>(T);
        tetra_fac<3, 4, 1>(T);
        tetra_fac<4, 5, 2>(T);
        tetra_fac<5, 0, 3>(T);
      }
      if (groups::group_IV_a) {
        tetra_fac<0, 2, 4>(T);
      }
      if (groups::group_IV_b) {
        tetra_fac<1, 3, 5>(T);
      }
      if (groups::group_V) {
        tetra_fac<0, 1, 6>(T);
        tetra_fac<1, 2, 6>(T);
        tetra_fac<2, 3, 6>(T);
        tetra_fac<3, 4, 6>(T);
        tetra_fac<4, 5, 6>(T);
        tetra_fac<5, 0, 6>(T);
      }
      if (groups::group_VI_a) {
        tetra_fac<0, 2, 6>(T);
        tetra_fac<2, 4, 6>(T);
        tetra_fac<4, 0, 6>(T);
      }
      if (groups::group_VI_b) {
        tetra_fac<1, 3, 6>(T);
        tetra_fac<3, 5, 6>(T);
        tetra_fac<5, 1, 6>(T);
      }
    }

    /**
     * Triangle updates:
     */
    for (octant = 0; octant < 8; ++octant) {
      inds = oct2inds[octant];
      if (groups::do_tri11_updates) {
        tri_fac<0, 2>(T);
        tri_fac<2, 4>(T);
        tri_fac<4, 0>(T);
      }
      if (groups::do_tri12_updates) {
        tri_fac<0, 1>(T);
        tri_fac<2, 1>(T);
        tri_fac<2, 3>(T);
        tri_fac<4, 3>(T);
        tri_fac<4, 5>(T);
        tri_fac<0, 5>(T);
      }
      if (groups::do_tri13_updates) {
        tri_fac<0, 6>(T);
        tri_fac<2, 6>(T);
        tri_fac<4, 6>(T);
      }
      if (groups::do_tri22_updates) {
        tri_fac<1, 3>(T);
        tri_fac<3, 5>(T);
        tri_fac<5, 1>(T);
      }
      if (groups::do_tri23_updates) {
        tri_fac<1, 6>(T);
        tri_fac<3, 6>(T);
        tri_fac<5, 6>(T);
      }
    }
  }
  else {
    /**
     * Tetrahedron updates:
     */
    for (octant = 0; octant < 8; ++octant) {
      inds = oct2inds[octant];
      if (groups::group_I) {
        tetra<1, 2, 3, P011, P010, P110>(T);
        tetra<3, 4, 5, P110, P100, P101>(T);
        tetra<5, 0, 1, P101, P001, P011>(T);
      }
      if (groups::group_II) {
        tetra<0, 1, 3, P001, P011, P110>(T);
        tetra<1, 2, 4, P011, P010, P100>(T);
        tetra<2, 3, 5, P010, P110, P101>(T);
        tetra<3, 4, 0, P110, P100, P001>(T);
        tetra<4, 5, 1, P100, P101, P011>(T);
        tetra<5, 0, 2, P101, P001, P010>(T);
      }
      if (groups::group_III) {
        tetra<0, 1, 4, P001, P011, P100>(T);
        tetra<1, 2, 5, P011, P010, P101>(T);
        tetra<2, 3, 0, P010, P110, P001>(T);
        tetra<3, 4, 1, P110, P100, P011>(T);
        tetra<4, 5, 2, P100, P101, P010>(T);
        tetra<5, 0, 3, P101, P001, P110>(T);
      }
      if (groups::group_IV_a) {
        tetra<0, 2, 4, P001, P010, P100>(T);
      }
      if (groups::group_IV_b) {
        tetra<1, 3, 5, P011, P110, P101>(T);
      }
      if (groups::group_V) {
        tetra<0, 1, 6, P001, P011, P111>(T);
        tetra<1, 2, 6, P011, P010, P111>(T);
        tetra<2, 3, 6, P010, P110, P111>(T);
        tetra<3, 4, 6, P110, P100, P111>(T);
        tetra<4, 5, 6, P100, P101, P111>(T);
        tetra<5, 0, 6, P101, P001, P111>(T);
      }
      if (groups::group_VI_a) {
        tetra<0, 2, 6, P001, P010, P111>(T);
        tetra<2, 4, 6, P010, P100, P111>(T);
        tetra<4, 0, 6, P100, P001, P111>(T);
      }
      if (groups::group_VI_b) {
        tetra<1, 3, 6, P011, P110, P111>(T);
        tetra<3, 5, 6, P110, P101, P111>(T);
        tetra<5, 1, 6, P101, P011, P111>(T);
      }
    }

    /**
     * Triangle updates:
     */
    for (octant = 0; octant < 8; ++octant) {
      inds = oct2inds[octant];
      if (groups::do_tri11_updates) {
        tri<0, 2, P001, P010>(T);
        tri<2, 4, P010, P100>(T);
        tri<4, 0, P100, P001>(T);
      }
      if (groups::do_tri12_updates) {
        tri<0, 1, P001, P011>(T);
        tri<2, 1, P010, P011>(T);
        tri<2, 3, P010, P110>(T);
        tri<4, 3, P100, P110>(T);
        tri<4, 5, P100, P101>(T);
        tri<0, 5, P001, P101>(T);
      }
      if (groups::do_tri13_updates) {
        tri<0, 6, P001, P111>(T);
        tri<2, 6, P010, P111>(T);
        tri<4, 6, P100, P111>(T);
      }
      if (groups::do_tri22_updates) {
        tri<1, 3, P011, P110>(T);
        tri<3, 5, P110, P101>(T);
        tri<5, 1, P101, P011>(T);
      }
      if (groups::do_tri23_updates) {
        tri<1, 6, P011, P111>(T);
        tri<3, 6, P110, P111>(T);
        tri<5, 6, P101, P111>(T);
      }
    }
  }

#if PRINT_UPDATES
  printf("olim3d::update_impl: T <- %g\n", T);
#endif
}

template <cost_func F, class node, int lp_norm, int d1, int d2>
void olim3d_hu<F, node, lp_norm, d1, d2>::update_crtp(double & T)
{
  using std::min;

  int i = n->get_i(), j = n->get_j(), k = n->get_k();
#if PRINT_UPDATES
  printf("olim3d_hu::update_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  double s_fac = INF(double);

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
  double p0[3], p1[3], p2[3], p_fac[3];

  // Depending on the compilation flag HU_USE_PARENT_NODE, either let
  // l0 be the index of the updating parent, or let l0 be the index of
  // the neighboring valid node with the minimum degree 0 update
  // value.
#if HU_USE_PARENT_NODE
  l0 = parent;
  p0[0] = __di(l0);
  p0[1] = __dj(l0);
  p0[2] = __dk(l0);
  T0 = updates::line<F, 3>()(//p0, VAL(l0), SPEED_ARGS(l0), h);
    p0, this->nb[l0]->get_value(), this->s_hat, this->s[l0], this->get_h());
#  if COLLECT_STATS
  node_stats.inc_line_updates(p0, 3);
#  endif
#else // HU_USE_PARENT_NODE
  for (int l = 0; l < 26; ++l) {
    if (nb[l]) {
      p0[0] = __di(l);
      p0[1] = __dj(l);
      p0[2] = __dk(l);
      // Tnew = this->template line<3>(p0, VAL(l), SPEED_ARGS(l), h);
      Tnew = updates::line<F, 3>()(
        p0, this->nb[l]->get_value(), this->s_hat, this->s[l], this->get_h());
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

  if (n->has_fac_parent()) {
    auto n_fac = static_cast<node *>(n->get_fac_parent());
    int i_fac = n_fac->get_i(), j_fac = n_fac->get_j(), k_fac = n_fac->get_k();
    s_fac = this->get_speed(i_fac, j_fac, k_fac);
    p_fac[0] = i_fac - i;
    p_fac[1] = j_fac - j;
    p_fac[2] = k_fac - k;
  }

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

      // TODO: skip triangle updates using KKT theory

      // Do the triangle update.
      auto const tmp = n->has_fac_parent() ?
        // this->template tri<3>(
        updates::tri<F, 3>()(
          // p0, p1, VAL(l0), VAL(l), SPEED_ARGS(l0, l), h, p_fac, s_fac)
          p0, p1, this->nb[l0]->get_value(), this->nb[l]->get_value(),
          this->s_hat, this->s[l0], this->s[l], this->get_h(), p_fac, s_fac) :
        updates::tri<F, 3>()(
          // p0, p1, VAL(l0), VAL(l), SPEED_ARGS(l0, l), h);
          p0, p1, this->nb[l0]->get_value(), this->nb[l]->get_value(),
          this->s_hat, this->s[l0], this->s[l], this->get_h());
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

        double mu[2], lam[2], df[2], d2f[3];
        int k;

        // TODO: Masha has some simplified way of computing the
        // Lagrange multipliers in the 3D Qpot paper---check this out
        auto const compute_lagrange_multipliers = [&] () {
          if (n->has_fac_parent()) {
            // F_fac_wkspc<F, 2> fac_wkspc; // TODO: enable
            assert(false);
          } else {
            F_wkspc<F, 2> w;
            set_args<F, 3>(
              w, p0, p1, p2,
              this->nb[l0]->get_value(),
              this->nb[l1]->get_value(),
              this->nb[l2]->get_value(),
              this->s_hat,
              this->s[l0],
              this->s[l1],
              this->s[l2],
              this->get_h());
            set_lambda<F, 3>(w, p0, p1, p2, lam);
            grad<2>(w, df);
            hess<2>(w, d2f);
            lagmults<2>(lam, df, d2f, mu, &k);
          }
        };

        // We get the first two checks for free using arglam
        assert(arglam[l1] != -1);
        lam[0] = arglam[l1];
        lam[1] = 0;
        compute_lagrange_multipliers();
        if (mu[0] < 0 || (k == 2 && mu[1] < 0)) {
          continue;
        }

        // TODO: not totally sure if this should be -1 or if it might
        // be -1 in some cases...
        if (arglam[l2] != -1) {
          lam[0] = 0;
          lam[1] = arglam[l2];
          compute_lagrange_multipliers();
          if (mu[0] < 0 || (k == 2 && mu[1] < 0)) {
            continue;
          }
        }

        // TODO: We're not doing the third triangle update right
        // now. This seems to work okay for now, but we can't do the
        // "exact solve" using the QR decomposition if we don't do the
        // third update. It's unclear how efficient this will be...

        // Finally, do the tetrahedron update.
        auto const tmp = n->has_fac_parent() ?
          updates::tetra<F, 3>()(
            //p0, p1, p2, VAL(l0), VAL(l1), VAL(l2), SPEED_ARGS(l0, l1, l2), h, 
            // p_fac, s_fac) :
            p0, p1, p2,
            this->nb[l0]->get_value(),
            this->nb[l1]->get_value(),
            this->nb[l2]->get_value(),
            this->s_hat, this->s[l0], this->s[l1], this->s[l2],
            this->get_h(), p_fac, s_fac) :
          updates::tetra<F, 3>()(
            // p0, p1, p2, VAL(l0), VAL(l1), VAL(l2), SPEED_ARGS(l0, l1, l2), h);
            p0, p1, p2,
            this->nb[l0]->get_value(),
            this->nb[l1]->get_value(),
            this->nb[l2]->get_value(),
            this->s_hat, this->s[l0], this->s[l1], this->s[l2],
            this->get_h());
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
          class tetra_updates, class groups, int num_neighbors>
olim3d_node_stats &
abstract_olim3d<
  node, line_updates, tri_updates, tetra_updates, groups,
  num_neighbors>::get_node_stats(int i, int j, int k)
{
  return _node_stats[this->get_height()*(this->get_width()*k + j) + i];
}

template <class node, class line_updates, class tri_updates,
          class tetra_updates, class groups, int num_neighbors>
olim3d_node_stats const &
abstract_olim3d<
  node, line_updates, tri_updates, tetra_updates, groups,
  num_neighbors>::get_node_stats(int i, int j, int k) const
{
  return _node_stats[this->get_height()*(this->get_width()*k + j) + i];
}

#endif // COLLECT_STATS

#endif // __OLIM3D_IMPL_HPP__
