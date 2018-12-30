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

template <class base, class node, int num_neighbors>
void abstract_olim3d<base, node, num_neighbors>::init()
{
  static_cast<base *>(this)->init_crtp();
}

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
      this->s[l] = this->get_speed(i + di<3>[l], j + dj<3>[l], k + dk<3>[l]);
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

#if COLLECT_STATS
  node_stats = this->get_node_stats(i, j, k);
  node_stats.inc_num_visits();
#endif

  // Do line update corresponding to parent node.
  {
    double Tnew = inf<double>;
    if (parent < 6) {
      line<1>(parent, Tnew);
    } else if (groups::num_neighbors > 6 && parent < 18) {
      line<2>(parent, Tnew);
    } else if (groups::num_neighbors > 18) {
      line<3>(parent, Tnew);
    }
    assert(!isinf(Tnew));
#if TRACK_PARENTS
    if (Tnew < T) {
      T = Tnew;
      n->set_parents({{
        &this->operator()(i + di<3>[parent], j + dj<3>[parent], k + dk<3>[parent]),
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
  reset_tri_skip_list();

  if (n->has_fac_parent()) {
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
}

template <cost_func F, class node, int lp_norm, int d1, int d2>
void olim3d_hu<F, node, lp_norm, d1, d2>::init_crtp()
{
  valid_d1 = new bool[26*26];
  valid_d2 = new bool[26*26];
  coplanar = new bool[26*26*26];

  static constexpr double tol = eps<double>;

  auto const is_valid = [&] (int l0, int l1, int d) -> bool {
    get_p(l0, p0);
    get_p(l1, p1);
    if (lp_norm == L1) {
      return dist1<3>(p0, p1) <= d + tol;
    } else if (lp_norm == L2) {
      return dist2sq<3>(p0, p1) <= d + tol;
    } else {
      return distmax<3>(p0, p1) <= d + tol;
    }
  };

  for (int l0 = 0; l0 < 26; ++l0) {
    for (int l1 = 0; l1 < 26; ++l1) {
      is_valid_d1(l0, l1) = is_valid(l0, l1, d1);
    }
  }

  for (int l0 = 0; l0 < 26; ++l0) {
    for (int l1 = 0; l1 < 26; ++l1) {
      is_valid_d2(l0, l1) = is_valid(l0, l1, d2);
    }
  }

  // Use the scalar triple product (dot(p x q, r)) to check if
  // three points are coplanar.
  {
    double p0_x_p1[3];
    for (int l0 = 0; l0 < 26; ++l0) {
      get_p(l0, p0);
      for (int l1 = 0; l1 < 26; ++l1) {
        get_p(l1, p1);
        p0_x_p1[0] = p0[1]*p1[2] - p0[2]*p1[1];
        p0_x_p1[1] = p0[2]*p1[0] - p0[0]*p1[2];
        p0_x_p1[2] = p0[0]*p1[1] - p0[1]*p1[0];
        for (int l2 = 0; l2 < 26; ++l2) {
          get_p(l2, p2);
          is_coplanar(l0, l1, l2) = fabs(
            p0_x_p1[0]*p2[0] + p0_x_p1[1]*p2[1] + p0_x_p1[2]*p2[2]) < 1e1*tol;
        }
      }
    }
  }
}

template <cost_func F, class node, int lp_norm, int d1, int d2>
void olim3d_hu<F, node, lp_norm, d1, d2>::update_crtp(double & T)
{
  using std::min;

  int i = n->get_i(), j = n->get_j(), k = n->get_k();
#if COLLECT_STATS
  auto & node_stats = this->get_node_stats(i, j, k);
  node_stats.inc_num_visits();
#endif

  /**
   * Tnew: temporary variable used for updating T, the update value
   * T0, T1, T2: separate minimum values for degree 0/1/2 updates
   * l0, l1: indices of p0 and p1 in 26 point neighborhood
   * p0, p1, p2: update node vectors
   */
  double Tnew, T0, T1 = inf<double>, T2 = inf<double>;
  int l0 = parent, l1 = -1;
  double s_fac;

  get_p(l0, p0);

  T0 = updates::line<F>()(
    get_p_norm(l0), this->nb[l0]->get_value(), this->s_hat, this->s[l0],
    this->get_h());
#if COLLECT_STATS
  node_stats.inc_line_updates(p0, 3);
#endif

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
  std::fill(arglam, arglam + 26, -1);

  // Find the minimal triangle update containing l0.
  for (int l = 0; l < 26; ++l) {
    if (l == l0 || !nb[l] || !is_valid_d1(l0, l)) {
      continue;
    }

    get_p(l, p1);

    // Do the triangle update.
    auto const tmp = n->has_fac_parent() ?
      updates::tri<F, 3>()(
        p0, p1, this->nb[l0]->get_value(), this->nb[l]->get_value(),
        this->s_hat, this->s[l0], this->s[l], this->get_h(), p_fac, s_fac) :
      updates::tri<F, 3>()(
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
  // There may only be a single valid neighbor, in which case we
  // should jump to the end of this function and skip all of the
  // tetrahedron updates.
  if (l1 == -1) {
    assert(isinf(T1));
    goto coda;
  } else {
    get_p(l1, p1);
  }

  // Do the tetrahedron updates such that p2 is sufficiently near p0
  // and p1.
  for (int l2 = 0; l2 < 26; ++l2) {
    if (l0 == l2 || l1 == l2 || !nb[l2] ||
        !is_valid_d2(l0, l2) || !is_valid_d2(l1, l2) ||
        is_coplanar(l0, l1, l2)) {
      continue;
    }

    get_p(l2, p2);

    updates::info<2> info;
    info.lambda[0] = arglam[l1];
    info.lambda[1] = 0;

    if (n->has_fac_parent()) {
      updates::tetra<F, 3>()(
        p0, p1, p2,
        this->nb[l0]->get_value(),
        this->nb[l1]->get_value(),
        this->nb[l2]->get_value(),
        this->s_hat, this->s[l0], this->s[l1], this->s[l2],
        this->get_h(), p_fac, s_fac, info);
    } else {
      updates::tetra<F, 3>()(
        p0, p1, p2,
        this->nb[l0]->get_value(),
        this->nb[l1]->get_value(),
        this->nb[l2]->get_value(),
        this->s_hat, this->s[l0], this->s[l1], this->s[l2],
        this->get_h(), info);
    }
#if COLLECT_STATS
    node_stats.inc_tetra_updates(p0, p1, p2, 3, tmp.is_degenerate(), true);
#endif
    Tnew = info.value;
    if (Tnew < T2) {
      T2 = Tnew;
    }
  }

  // Set T to be the minimum of T0, T1, and T2.
coda:
  T = min(T0, min(T1, T2));
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
