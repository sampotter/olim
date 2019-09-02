#pragma once

#include <src/config.hpp>

#include "offsets.hpp"

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

template <cost_func F, class groups, ordering ord>
void eikonal::olim3d_bv<F, groups, ord>::update_impl(
  int lin_hat, int const * nb, int parent, double & U)
{
  using std::min;

  assert(parent <= groups::num_nb);

  // Do line update corresponding to parent node.
  //
  // TODO: we can combine the one-point update done by the hu and bv
  // implementations and save ourselves a few lines of code.
  {
    double U_ = inf<double>;
    if (parent < 6) {
      line<1>(lin_hat, nb, parent, U_);
    } else if (groups::num_nb > 6 && parent < 18) {
      line<2>(lin_hat, nb, parent, U_);
    } else if (groups::num_nb > 18) {
      line<3>(lin_hat, nb, parent, U_);
    }
    assert(!isinf(U_));
    U = min(U, U_);
    // TODO: collect stats
  }

  /**
   * This array is used to keep track of which triangle updates should
   * be skipped.
   *
   * TODO: we could use a bitvector here instead.
   */
  reset_tri_skip_list();

  if (this->is_factored(lin_hat)) {
    /**
     * Tetrahedron updates:
     */
    // TODO: decide which octants to do based on the parent index
    // (only need to look at 1, 2, or 4 octants, never more than that)
    for (octant = 0; octant < 8; ++octant) {
      inds = oct2inds[octant];
      if (groups::group_I) {
        tetra_fac<1, 2, 3>(lin_hat, nb, parent, U);
        tetra_fac<3, 4, 5>(lin_hat, nb, parent, U);
        tetra_fac<5, 0, 1>(lin_hat, nb, parent, U);
      }
      if (groups::group_II) {
        tetra_fac<0, 1, 3>(lin_hat, nb, parent, U);
        tetra_fac<1, 2, 4>(lin_hat, nb, parent, U);
        tetra_fac<2, 3, 5>(lin_hat, nb, parent, U);
        tetra_fac<3, 4, 0>(lin_hat, nb, parent, U);
        tetra_fac<4, 5, 1>(lin_hat, nb, parent, U);
        tetra_fac<5, 0, 2>(lin_hat, nb, parent, U);
      }
      if (groups::group_III) {
        tetra_fac<0, 1, 4>(lin_hat, nb, parent, U);
        tetra_fac<1, 2, 5>(lin_hat, nb, parent, U);
        tetra_fac<2, 3, 0>(lin_hat, nb, parent, U);
        tetra_fac<3, 4, 1>(lin_hat, nb, parent, U);
        tetra_fac<4, 5, 2>(lin_hat, nb, parent, U);
        tetra_fac<5, 0, 3>(lin_hat, nb, parent, U);
      }
      if (groups::group_IV_a) {
        tetra_fac<0, 2, 4>(lin_hat, nb, parent, U);
      }
      if (groups::group_IV_b) {
        tetra_fac<1, 3, 5>(lin_hat, nb, parent, U);
      }
      if (groups::group_V) {
        tetra_fac<0, 1, 6>(lin_hat, nb, parent, U);
        tetra_fac<1, 2, 6>(lin_hat, nb, parent, U);
        tetra_fac<2, 3, 6>(lin_hat, nb, parent, U);
        tetra_fac<3, 4, 6>(lin_hat, nb, parent, U);
        tetra_fac<4, 5, 6>(lin_hat, nb, parent, U);
        tetra_fac<5, 0, 6>(lin_hat, nb, parent, U);
      }
      if (groups::group_VI_a) {
        tetra_fac<0, 2, 6>(lin_hat, nb, parent, U);
        tetra_fac<2, 4, 6>(lin_hat, nb, parent, U);
        tetra_fac<4, 0, 6>(lin_hat, nb, parent, U);
      }
      if (groups::group_VI_b) {
        tetra_fac<1, 3, 6>(lin_hat, nb, parent, U);
        tetra_fac<3, 5, 6>(lin_hat, nb, parent, U);
        tetra_fac<5, 1, 6>(lin_hat, nb, parent, U);
      }
    }

    /**
     * Triangle updates:
     */
    for (octant = 0; octant < 8; ++octant) {
      inds = oct2inds[octant];
      if (groups::do_tri11_updates) {
        tri_fac<0, 2>(lin_hat, nb, parent, U);
        tri_fac<2, 4>(lin_hat, nb, parent, U);
        tri_fac<4, 0>(lin_hat, nb, parent, U);
      }
      if (groups::do_tri12_updates) {
        tri_fac<0, 1>(lin_hat, nb, parent, U);
        tri_fac<2, 1>(lin_hat, nb, parent, U);
        tri_fac<2, 3>(lin_hat, nb, parent, U);
        tri_fac<4, 3>(lin_hat, nb, parent, U);
        tri_fac<4, 5>(lin_hat, nb, parent, U);
        tri_fac<0, 5>(lin_hat, nb, parent, U);
      }
      if (groups::do_tri13_updates) {
        tri_fac<0, 6>(lin_hat, nb, parent, U);
        tri_fac<2, 6>(lin_hat, nb, parent, U);
        tri_fac<4, 6>(lin_hat, nb, parent, U);
      }
      if (groups::do_tri22_updates) {
        tri_fac<1, 3>(lin_hat, nb, parent, U);
        tri_fac<3, 5>(lin_hat, nb, parent, U);
        tri_fac<5, 1>(lin_hat, nb, parent, U);
      }
      if (groups::do_tri23_updates) {
        tri_fac<1, 6>(lin_hat, nb, parent, U);
        tri_fac<3, 6>(lin_hat, nb, parent, U);
        tri_fac<5, 6>(lin_hat, nb, parent, U);
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
        tetra<1, 2, 3, 0b011, 0b010, 0b110>(lin_hat, nb, parent, U);
        tetra<3, 4, 5, 0b110, 0b100, 0b101>(lin_hat, nb, parent, U);
        tetra<5, 0, 1, 0b101, 0b001, 0b011>(lin_hat, nb, parent, U);
      }
      if (groups::group_II) {
        tetra<0, 1, 3, 0b001, 0b011, 0b110>(lin_hat, nb, parent, U);
        tetra<1, 2, 4, 0b011, 0b010, 0b100>(lin_hat, nb, parent, U);
        tetra<2, 3, 5, 0b010, 0b110, 0b101>(lin_hat, nb, parent, U);
        tetra<3, 4, 0, 0b110, 0b100, 0b001>(lin_hat, nb, parent, U);
        tetra<4, 5, 1, 0b100, 0b101, 0b011>(lin_hat, nb, parent, U);
        tetra<5, 0, 2, 0b101, 0b001, 0b010>(lin_hat, nb, parent, U);
      }
      if (groups::group_III) {
        tetra<0, 1, 4, 0b001, 0b011, 0b100>(lin_hat, nb, parent, U);
        tetra<1, 2, 5, 0b011, 0b010, 0b101>(lin_hat, nb, parent, U);
        tetra<2, 3, 0, 0b010, 0b110, 0b001>(lin_hat, nb, parent, U);
        tetra<3, 4, 1, 0b110, 0b100, 0b011>(lin_hat, nb, parent, U);
        tetra<4, 5, 2, 0b100, 0b101, 0b010>(lin_hat, nb, parent, U);
        tetra<5, 0, 3, 0b101, 0b001, 0b110>(lin_hat, nb, parent, U);
      }
      if (groups::group_IV_a) {
        tetra<0, 2, 4, 0b001, 0b010, 0b100>(lin_hat, nb, parent, U);
      }
      if (groups::group_IV_b) {
        tetra<1, 3, 5, 0b011, 0b110, 0b101>(lin_hat, nb, parent, U);
      }
      if (groups::group_V) {
        tetra<0, 1, 6, 0b001, 0b011, 0b111>(lin_hat, nb, parent, U);
        tetra<1, 2, 6, 0b011, 0b010, 0b111>(lin_hat, nb, parent, U);
        tetra<2, 3, 6, 0b010, 0b110, 0b111>(lin_hat, nb, parent, U);
        tetra<3, 4, 6, 0b110, 0b100, 0b111>(lin_hat, nb, parent, U);
        tetra<4, 5, 6, 0b100, 0b101, 0b111>(lin_hat, nb, parent, U);
        tetra<5, 0, 6, 0b101, 0b001, 0b111>(lin_hat, nb, parent, U);
      }
      if (groups::group_VI_a) {
        tetra<0, 2, 6, 0b001, 0b010, 0b111>(lin_hat, nb, parent, U);
        tetra<2, 4, 6, 0b010, 0b100, 0b111>(lin_hat, nb, parent, U);
        tetra<4, 0, 6, 0b100, 0b001, 0b111>(lin_hat, nb, parent, U);
      }
      if (groups::group_VI_b) {
        tetra<1, 3, 6, 0b011, 0b110, 0b111>(lin_hat, nb, parent, U);
        tetra<3, 5, 6, 0b110, 0b101, 0b111>(lin_hat, nb, parent, U);
        tetra<5, 1, 6, 0b101, 0b011, 0b111>(lin_hat, nb, parent, U);
      }
    }

    /**
     * Triangle updates:
     */
    for (octant = 0; octant < 8; ++octant) {
      inds = oct2inds[octant];
      if (groups::do_tri11_updates) {
        tri<0, 2, 0b001, 0b010>(lin_hat, nb, parent, U);
        tri<2, 4, 0b010, 0b100>(lin_hat, nb, parent, U);
        tri<4, 0, 0b100, 0b001>(lin_hat, nb, parent, U);
      }
      if (groups::do_tri12_updates) {
        tri<0, 1, 0b001, 0b011>(lin_hat, nb, parent, U);
        tri<2, 1, 0b010, 0b011>(lin_hat, nb, parent, U);
        tri<2, 3, 0b010, 0b110>(lin_hat, nb, parent, U);
        tri<4, 3, 0b100, 0b110>(lin_hat, nb, parent, U);
        tri<4, 5, 0b100, 0b101>(lin_hat, nb, parent, U);
        tri<0, 5, 0b001, 0b101>(lin_hat, nb, parent, U);
      }
      if (groups::do_tri13_updates) {
        tri<0, 6, 0b001, 0b111>(lin_hat, nb, parent, U);
        tri<2, 6, 0b010, 0b111>(lin_hat, nb, parent, U);
        tri<4, 6, 0b100, 0b111>(lin_hat, nb, parent, U);
      }
      if (groups::do_tri22_updates) {
        tri<1, 3, 0b011, 0b110>(lin_hat, nb, parent, U);
        tri<3, 5, 0b110, 0b101>(lin_hat, nb, parent, U);
        tri<5, 1, 0b101, 0b011>(lin_hat, nb, parent, U);
      }
      if (groups::do_tri23_updates) {
        tri<1, 6, 0b011, 0b111>(lin_hat, nb, parent, U);
        tri<3, 6, 0b110, 0b111>(lin_hat, nb, parent, U);
        tri<5, 6, 0b101, 0b111>(lin_hat, nb, parent, U);
      }
    }
  }
}

template <cost_func F, int lp_norm, int d1, int d2, ordering ord>
void eikonal::olim3d_hu<F, lp_norm, d1, d2, ord>::init()
{
  // TODO: only allocate once
  valid_d1 = new bool[26*26];
  valid_d2 = new bool[26*26];
  coplanar = new bool[26*26*26];
  geom_wkspcs = new geom_wkspc<2>[26*26*26];
  qr_wkspcs = new qr_wkspc<3, 2>[26*26*26];

  static constexpr double tol = eps<double>;

  auto const is_valid = [&] (vec2<int> inds, int d) -> bool {
    p0 = get_p<3>(inds[0]);
    p1 = get_p<3>(inds[1]);
    if (lp_norm == L1) {
      return dist1(p0, p1) <= d + tol;
    } else if (lp_norm == L2) {
      return dist2sq(p0, p1) <= d + tol;
    } else {
      return disti(p0, p1) <= d + tol;
    }
  };

  for (auto inds: range<2> {{26, 26}}) {
    is_valid_d1(inds) = is_valid(inds, d1);
  }

  for (auto inds: range<2> {{26, 26}}) {
    is_valid_d2(inds) = is_valid(inds, d2);
  }

  // Use the scalar triple product (dot(p x q, r)) to check if
  // three points are coplanar.
  for (auto inds: range<3> {{26, 26, 26}}) {
    p0 = get_p<3>(inds[0]);
    p1 = get_p<3>(inds[1]);
    p2 = get_p<3>(inds[2]);
    is_coplanar(inds) = fabs((p0^p1)*p2) < 1e1*tol;
    get_geom_wkspc(inds).template init<3>(p0, p1, p2);
    get_qr_wkspc(inds).init(p0, p1, p2);
  }
}

template <cost_func F, int lp_norm, int d1, int d2, ordering ord>
void eikonal::olim3d_hu<F, lp_norm, d1, d2, ord>::update_impl(
  int lin_hat, int const * nb, int parent, double & U)
{
  using std::min;

  /**
   * U_: temporary variable used for updating U, the update value
   * U0, U1, U2: separate minimum values for degree 0/1/2 updates
   * l0, l1: indices of p0 and p1 in 26 point neighborhood
   * p0, p1, p2: update node vectors
   */
  double U_, U0, U1 = inf<double>, U2 = inf<double>;
  int l0 = parent, l1 = -1;
  double s_fac = inf<double>; // Silence warning

  p0 = get_p<3>(l0);

  // TODO: see comment above about one-point updates
  U0 = eikonal::line<F>()(
    get_p_norm(l0), this->_U[nb[l0]], this->_s[lin_hat], this->_s[nb[l0]],
    this->get_h());

  if (this->is_factored(lin_hat)) {
    auto src = this->_fac_srcs[lin_hat];
    s_fac = src->s;
    p_fac = vec3<double> {src->x} - this->to_vector_index(lin_hat);
  }

  // Create a cache for the minimizing lambdas to use for skipping
  // tetrahedron updates. Initialize it to -1 so that we can tell
  // which triangle updates have been computed.
  double arglam[26];
  std::fill(arglam, arglam + 26, -1);

  // Find the minimal triangle update containing l0.
  for (int l = 0; l < 26; ++l) {
    if (l == l0 || nb[l] == -1 || !is_valid_d1({l0, l})) {
      continue;
    }

    p1 = get_p<3>(l);

    // Do the triangle update.
    auto const tmp = this->is_factored(lin_hat) ?
      eikonal::tri<F, 3>()(
        p0, p1, this->_U[nb[l0]], this->_U[nb[l]],
        this->_s[lin_hat], this->_s[nb[l0]], this->_s[nb[l]], this->get_h(), p_fac, s_fac) :
      eikonal::tri<F, 3>()(
        p0, p1, this->_U[nb[l0]], this->_U[nb[l]],
        this->_s[lin_hat], this->_s[nb[l0]], this->_s[nb[l]], this->get_h());
    U_ = tmp.value;
    if (U_ < U1) {
      U1 = U_;
      l1 = l;
    }
    arglam[l] = tmp.lambda[0];
  }
  // There may only be a single valid neighbor, in which case we
  // should jump to the end of this function and skip all of the
  // tetrahedron updates.
  if (l1 == -1) {
    assert(isinf(U1));
    goto coda;
  } else {
    p1 = get_p<3>(l1);
  }

  // Do the tetrahedron updates such that p2 is sufficiently near p0
  // and p1.
  for (int l2 = 0; l2 < 26; ++l2) {
    if (l0 == l2 || l1 == l2 || nb[l2] == -1 ||
        !is_valid_d2({l0, l2}) || !is_valid_d2({l1, l2}) ||
        is_coplanar({l0, l1, l2})) {
      continue;
    }

    p2 = get_p<3>(l2);

    update_info<2> info;
    info.lambda[0] = arglam[l1];
    info.lambda[1] = 0;

    double u0 = this->_U[nb[l0]], u1 = this->_U[nb[l1]],
      u2 = this->_U[nb[l2]], s = this->_s[lin_hat], s0 = this->_s[nb[l0]],
      s1 = this->_s[nb[l1]], s2 = this->_s[nb[l2]], h = this->get_h();

    if (this->is_factored(lin_hat)) {
      geom_fac_wkspc<2> g;
      g.init<3>(p0, p1, p2, p_fac);
      F_fac_wkspc<F, 2> w;
      set_args<F>(w, g, u0, u1, u2, s, s0, s1, s2, h, s_fac);
      cost_functor_fac<F, 3, 2> func {w, g};
      eikonal::tetra<F, 3>()(func, info);
      if (F == MP0) {
        eval_mp1_fix(func.w, s, s0, s1, s2, h, info.lambda, info.value);
      }
    }
    else {
      F_wkspc<F, 2> w;
      set_args<F>(w, u0, u1, u2, s, s0, s1, s2, h);
      int lin = to_nb_linear_index(vec3<int> {l0, l1, l2});
      cost_functor<F, 3, 2> func {w, geom_wkspcs[lin]};
      func.qr = &qr_wkspcs[lin];
      eikonal::tetra<F, 3>()(func, info);
      if (F == MP0 && info.inbounds()) {
        func.set_lambda(info.lambda);
        eval_mp1_fix(func.w, s, s0, s1, s2, h, info.lambda, info.value);
      }
    }

    U_ = info.value;
    if (U_ < U2) {
      U2 = U_;
    }
  }

coda:
  U = min(U0, min(U1, U2));
}
