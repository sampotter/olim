#ifndef __OLIM_TEST_COMMON_HPP__
#define __OLIM_TEST_COMMON_HPP__

#include <gtest/gtest.h>

#include <cassert>
#include <cmath>
#include <type_traits>
#include <vector>

#include "common.hpp"
#include "range.hpp"
#include "slow.hpp"
#include "vec.hpp"

testing::AssertionResult
eq_w_rel_tol(double x, double y, double tol = 1e-13) {
  if (x == 0 && y == 0) {
    return testing::AssertionSuccess();
  } else {
    double denom = fmax(fabs(x), fabs(y));
    double rel_err = fabs(x - y)/denom;
    if (rel_err > tol) {
      return testing::AssertionFailure()
        << "|" << x << " - " << y << "|/|" << denom << "| > " << tol;
    } else {
      return testing::AssertionSuccess();
    }
  }
}

template <class olim>
testing::AssertionResult
trivial_case_works() {
  olim m {{1, 1}, 1};
  m.add_boundary_node({0, 0});
  m.run();
  return eq_w_rel_tol(m.get_U({0, 0}), 0.0);
}

template <class olim>
testing::AssertionResult
adjacent_update_works() {
  olim m {{2, 1}, 0.5};
  m.add_boundary_node({0, 0});
  m.run();
  return eq_w_rel_tol(m.get_U({1, 0}), 0.5);
}

template <class olim>
testing::AssertionResult
correct_corners_in_limit(int n, double tol) {
  double h = 2./(n - 1);
  olim m {{n, n}, h, s0<2>, {1., 1.}};
  m.add_boundary_node({n/2, n/2});
  m.run();

  auto res = eq_w_rel_tol(m.get_U({0, n - 1}), sqrt2, tol);
  if (!res) return res;
  
  res = eq_w_rel_tol(m.get_U({0, n - 1}), sqrt2, tol);
  if (!res) return res;

  res = eq_w_rel_tol(m.get_U({n - 1, 0}), sqrt2, tol);
  if (!res) return res;

  res = eq_w_rel_tol(m.get_U({n - 1, n - 1}), sqrt2, tol);
  if (!res) return res;

  return testing::AssertionSuccess();
}

template <class olim>
testing::AssertionResult
quadrants_are_correct(
  double diag_value, std::enable_if_t<olim::ndim == 2> * = 0)
{
  int n = 2;
  double h = 1;

  double x0[4] = {0.0, 0.0, 1.0, 1.0};
  double y0[4] = {0.0, 1.0, 0.0, 1.0};

  for (int i = 0, i0 = 0, j0 = 0; i < 4; i0 = ++i/2, j0 = i0 % 2) {
    olim m {{n, n}, h, s0<2>, {x0[i], y0[i]}};
    m.add_boundary_node({i0, j0});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({i0, j0}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({(i0 + 1) % 2, j0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({i0, (j0 + 1) % 2}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({(i0 + 1) % 2, (j0 + 1) % 2}), diag_value);
    if (!res) return res;
  }

  return testing::AssertionSuccess();
}

template <class olim3d_t>
testing::AssertionResult
quadrants_are_correct(
  double diag_value, std::enable_if_t<olim3d_t::ndim == 3> * = 0)
{
  int n = 2;
  double h = 1;

  /*
   * Tests for quadrants in the x-y plane:
   */
  {
    olim3d_t m {{n, n, 1}, h, s0<3>, {0, 0, 0}};
    m.add_boundary_node({0, 0, 0});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 1, 0}), diag_value);
    if (!res) return res;
  }
  {
    olim3d_t m {{n, n, 1}, h, s0<3>, {1, 0, 0}};
    m.add_boundary_node({0, 1, 0});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 0}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 0}), diag_value);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 1, 0}), 1.0);
    if (!res) return res;
  }
  {
    olim3d_t m {{n, n, 1}, h, s0<3>, {0, 1, 0}};
    m.add_boundary_node({1, 0, 0});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 0}), diag_value);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 0}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 1, 0}), 1.0);
    if (!res) return res;
  }
  {
    olim3d_t m {{n, n, 1}, h, s0<3>, {1, 1, 0}};
    m.add_boundary_node({1, 1, 0});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), diag_value);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 1, 0}), 0.0);
    if (!res) return res;
  }

  /**
   * Tests for quadrants in the x-z plane:
   */
  {
    olim3d_t m {{n, 1, n}, h, s0<3>, {0, 0, 0}};
    m.add_boundary_node({0, 0, 0});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 0, 1}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 1}), diag_value);
    if (!res) return res;
  }
  {
    olim3d_t m {{n, 1, n}, h, s0<3>, {0, 1, 0}};
    m.add_boundary_node({1, 0, 0});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 0, 1}), diag_value);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 0}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 1}), 1.0);
    if (!res) return res;
  }
  {
    olim3d_t m {{n, 1, n}, h, s0<3>, {0, 0, 1}};
    m.add_boundary_node({0, 0, 1});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 0, 1}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 0}), diag_value);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 1}), 1.0);
    if (!res) return res;
  }
  {
    olim3d_t m {{n, 1, n}, h, s0<3>, {0, 1, 1}};
    m.add_boundary_node({1, 0, 1});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), diag_value);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 0, 1}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({1, 0, 1}), 0.0);
    if (!res) return res;
  }

  /**
   * Tests for quadrants in the y-z plane:
   */
  {
    olim3d_t m {{1, n, n}, h, s0<3>, {0, 0, 0}};
    m.add_boundary_node({0, 0, 0});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 0, 1}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 1}), diag_value);
    if (!res) return res;
  }
  {
    olim3d_t m {{1, n, n}, h, s0<3>, {1, 0, 0}};
    m.add_boundary_node({0, 1, 0});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 0, 1}), diag_value);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 0}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 1}), 1.0);
    if (!res) return res;
  }
  {
    olim3d_t m {{1, n, n}, h, s0<3>, {0, 0, 1}};
    m.add_boundary_node({0, 0, 1});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 0, 1}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 0}), diag_value);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 1}), 1.0);
    if (!res) return res;
  }
  {
    olim3d_t m {{1, n, n}, h, s0<3>, {1, 0, 1}};
    m.add_boundary_node({0, 1, 1});
    m.run();

    auto res = eq_w_rel_tol(m.get_U({0, 0, 0}), diag_value);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 0, 1}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 0}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({0, 1, 1}), 0.0);
    if (!res) return res;
  }

  return testing::AssertionSuccess();
}

template <class olim3d_t>
testing::AssertionResult
octants_are_correct(double diag2val, double diag3val) {
  vec3<int> dims {2, 2, 2};

  for (auto inds: range<3> {dims}) {
    olim3d_t m {dims, 1., s0<3>, vec3<double> {inds}};
    m.add_boundary_node(inds);
    m.run();

    int i = inds[0], j = inds[1], k = inds[2];

    auto res = eq_w_rel_tol(m.get_U({i, j, k}), 0.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({(i + 1) % 2, j, k}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({i, (j + 1) % 2, k}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({i, j, (k + 1) % 2}), 1.0);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({(i + 1) % 2, (j + 1) % 2, k}), diag2val);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({(i + 1) % 2, j, (k + 1) % 2}), diag2val);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({i, (j + 1) % 2, (k + 1) % 2}), diag2val);
    if (!res) return res;

    res = eq_w_rel_tol(m.get_U({(i + 1) % 2, (j + 1) % 2, (k + 1) % 2}), diag3val);
    if (!res) return res;
  }

  return testing::AssertionSuccess();
}

template <class olim, class olim3d_t>
testing::AssertionResult
planes_are_correct(slow<2> s = s0<2>, slow<3> s3d = s0<3>, int n = 51)
{
  double h = 2.0/(n - 1);
  
  olim m {{n, n}, h, s, {1, 1}};
  m.add_boundary_node({n/2, n/2});
  m.run();
  
  olim3d_t m3d {{n, n, n}, h, s3d, {1, 1, 1}};
  m3d.add_boundary_node({n/2, n/2, n/2});
  m3d.run();

  auto msg = [] (testing::AssertionResult & res, int i, int j, int k)
    -> testing::AssertionResult &
  {
    return res << ", (i = " << i << ", j = " << j << ", k = " << k << ")";
  };

  for (auto inds: range<2> {{n, n}}) {
    double U = m.get_U(inds);
    int i = inds[0], j = inds[1];

    auto res = eq_w_rel_tol(U, m3d.get_U({i, j, n/2}));
    if (!res) return msg(res, i, j, n/2);

    res = eq_w_rel_tol(U, m3d.get_U({i, n/2, j}));
    if (!res) return msg(res, i, n/2, j);

    res = eq_w_rel_tol(U, m3d.get_U({n/2, i, j}));
    if (!res) return msg(res, n/2, i, j);
  }

  return testing::AssertionSuccess();
}

template <class olim>
testing::AssertionResult
result_is_symmetric(slow<2> s = s0, int n = 51, double tol = 1e-13) {
  double h = 2.0/(n - 1);
  olim m {{n, n}, h, s, {1, 1}};
  m.add_boundary_node({n/2, n/2});
  m.run();

  auto msg = [] (testing::AssertionResult & res, int i, int j)
    -> testing::AssertionResult &
  {
    return res << ", (i = " << i << ", j = " << j << ")";
  };

  for (int i = 0; i < n; ++i) {
    for (int j = 0, j_ = n - 1; j < n; ++j, --j_) {
      auto res = eq_w_rel_tol(m.get_U({i, j}), m.get_U({i, j_}), tol);
      if (!res) return msg(res, i, j);
    }
  }

  for (int i = 0, i_ = n - 1; i < n; ++i, --i_) {
    for (int j = 0; j < n; ++j) {
      auto res = eq_w_rel_tol(m.get_U({i, j}), m.get_U({i_, j}), tol);
      if (!res) return msg(res, i, j);
    }
  }

  return testing::AssertionSuccess();
}

template <class olim3d_t>
testing::AssertionResult
result_is_symmetric(slow<3> s = s0, int n = 21, double tol = 1e-13) {
  double h = 2.0/(n - 1);
  olim3d_t m {{n, n, n}, h, s, {1, 1, 1}};
  m.add_boundary_node({n/2, n/2, n/2});
  m.run();

  auto msg = [] (testing::AssertionResult & res, int i, int j, int k)
    -> testing::AssertionResult &
  {
    return res << ", (i = " << i << ", j = " << j << ", k = " << k << ")";
  };

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0, k_ = n - 1; k < n; ++k, --k_) {
        auto res = eq_w_rel_tol(
          m.get_U({i, j, k}), m.get_U({i, j, k_}), tol);
        if (!res) return msg(res, i, j, k);
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      for (int j = 0, j_ = n - 1; j < n; ++j, --j_) {
        auto res = eq_w_rel_tol(
          m.get_U({i, j, k}), m.get_U({i, j_, k}), tol);
        if (!res) return msg(res, i, j, k);
      }
    }
  }

  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) {
      for (int i = 0, i_ = n - 1; i < n; ++i, --i_) {
        auto res = eq_w_rel_tol(
          m.get_U({i, j, k}), m.get_U({i_, j, k}), tol);
        if (!res) return msg(res, i, j, k);
      }
    }
  }

  return testing::AssertionSuccess();
}

template <class olim3d_t>
testing::AssertionResult
two_by_two_by_three_cells_are_correct() {
  int dims[3][3] = {{3, 2, 2}, {2, 3, 2}, {2, 2, 3}};

  olim3d_t m_gt {{3, 2, 2}, 1, s0<3>, {0, 0, 0}};
  m_gt.add_boundary_node({0, 0, 0});
  m_gt.run();

  for (int dim = 0, imax, jmax, kmax; dim < 3; ++dim) {
    imax = dims[dim][0];
    jmax = dims[dim][1];
    kmax = dims[dim][2];

    for (int corner = 0; corner < 8; ++corner) {
      int i0 = (imax - 1)*(corner & 1);
      int j0 = (jmax - 1)*((corner & 2) >> 1);
      int k0 = (kmax - 1)*((corner & 4) >> 2);

      double x0 = j0, y0 = i0, z0 = k0;
      olim3d_t m {{imax, jmax, kmax}, 1, s0<3>, {x0, y0, z0}};
      m.add_boundary_node({i0, j0, k0});
      m.run();

      int da = i0 == 0 ? 1 : -1;
      int db = j0 == 0 ? 1 : -1;
      int dc = k0 == 0 ? 1 : -1;
      for (int i = 0, a = i0; i < imax; ++i, a += da) {
        for (int j = 0, b = j0; j < jmax; ++j, b += db) {
          for (int k = 0, c = k0, a_, b_, c_; k < kmax; ++k, c += dc) {
            a_ = imax == 3 ? a : jmax == 3 ? b : c;
            b_ = imax == 3 ? b : jmax == 3 ? c : a;
            c_ = imax == 3 ? c : jmax == 3 ? a : b;
            auto res = eq_w_rel_tol(m.get_U({i, j, k}), m_gt.get_U({a_, b_, c_}));
            if (!res) return res;
          }
        }
      }
    }
  }

  return testing::AssertionSuccess();
}

template <class olim3d_t>
testing::AssertionResult
plane_boundaries_are_correct() {
  int n = 2;
  double h = 1;
  olim3d_t m {{n, n, n}, h, s0<3>, {0, 0, 0}};

  vec3<int> inds[4];
  double Us[4];

  // typename olim3d_t::node_type nodes[4];
  int k = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      // nodes[k++] = typename olim3d_t::node_type {i, j, 0};
      inds[k] = {i, j, 0};
      Us[k++] = 0;
    }
  }
  m.add_boundary_nodes(inds, Us, 4);

  m.run();
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      auto res = eq_w_rel_tol(m.get_U({i, j, 1}), 1.0);
      if (!res) return res;
    }
  }

  return testing::AssertionSuccess();
}

template <class olim1, class olim2>
testing::AssertionResult
olims_agree(slow<2> s = s0, int n = 51) {
  double h = 2.0/(n - 1);
  int i0 = n/2;

  olim1 m1 {{n, n}, h, s, {1, 1}};
  m1.add_boundary_node({i0, i0});
  m1.run();

  olim2 m2 {{n, n}, h, s, {1, 1}};
  m2.add_boundary_node({i0, i0});
  m2.run();

  for (auto inds: range<2> {{n, n}}) {
    auto res = eq_w_rel_tol(m1.get_U(inds), m2.get_U(inds));
    if (!res) return res;
  }

  return testing::AssertionSuccess();
}

template <class olim3d_t, class other_olim3d_t>
void agrees_with_other_olim3d_t(int n = 21) {
  double h = 2.0/(n - 1);
  int i0 = n/2;

  olim3d_t m1 {n, n, n, h, s0<3>, 1, 1, 1};
  m1.add_boundary_node(i0, i0, i0);
  m1.run();

  other_olim3d_t m2 {n, n, n, h, s0<3>, 1, 1, 1};
  m2.add_boundary_node(i0, i0, i0);
  m2.run();

  for (auto inds: range<3> {{n, n, n}}) {
    auto res = eq_w_rel_tol(m1.get_U(inds), m2.get_U(inds));
    if (!res) return res;
  }
}

template <class olim, class olim3d_t>
testing::AssertionResult
planes_are_correct_nonsymmetric(
  double (* s3d)(double, double, double),
  double (* fxy)(double, double),
  double (* fyz)(double, double),
  double (* fxz)(double, double),
  int n1 = 15,
  int n2 = 20,
  int n3 = 25,
  int i1 = 5,
  int i2 = 10,
  int i3 = 15,
  double h = 1,
  double diam = 2)
{
  double x0 = (n2 - 1)/diam, y0 = (n1 - 1)/diam, z0 = (n3 - 1)/diam;

  olim3d_t m3d {n1, n2, n3, h, s3d, x0, y0, z0};
  m3d.add_boundary_node(i1, i2, i3);
  m3d.run();

  for (int i = 0; i < n1; ++i) {
    double y = i*h - y0;
    for (int j = 0; j < n2; ++j) {
      double x = j*h - x0;
      double u = fxy(x, y);
      auto res = eq_w_rel_tol(m3d.get_U(i, j, 0), u);
      if (!res) return res;
    }
  }

  for (int i = 0; i < n1; ++i) {
    double y = i*h - y0;
    for (int k = 0; k < n3; ++k) {
      double z = k*h - z0;
      double u = fyz(y, z);
      auto res = eq_w_rel_tol(m3d.get_U(i, 0, k), u);
      if (!res) return res;
    }
  }

  for (int j = 0; j < n2; ++j) {
    double x = j*h - x0;
    for (int k = 0; k < n3; ++k) {
      double z = k*h - z0;
      double u = fxz(x, z);
      auto res = eq_w_rel_tol(m3d.get_U(0, j, k), u);
      if (!res) return res;
    }
  }

  return testing::AssertionSuccess();
}

template <class olim, class olim3d_t>
testing::AssertionResult
planes_agree_nonsymmetric(
  double (* s3d)(double, double, double),
  double (* sxy)(double, double),
  double (* syz)(double, double),
  double (* sxz)(double, double),
  int n = 21)
{
  double h = 2.0/(n - 1);
  int i0 = n/2;

  olim3d_t m3d {n, n, n, h, s3d, 1, 1, 1};
  m3d.add_boundary_node(i0, i0, i0);
  m3d.run();

  olim mxy {n, n, h, sxy, 1, 1};
  mxy.add_boundary_node(i0, i0);
  mxy.run();

  olim myz {n, n, h, syz, 1, 1};
  myz.add_boundary_node(i0, i0);
  myz.run();

  olim mxz {n, n, h, sxz, 1, 1};
  mxz.add_boundary_node(i0, i0);
  mxz.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      auto res = eq_w_rel_tol(m3d.get_U(i, j, 0), mxy.get_U(i, j));
      if (!res) return res;
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      auto res = eq_w_rel_tol(m3d.get_U(i, 0, k), myz.get_U(i, k));
      if (!res) return res;
    }
  }

  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) {
      auto res = eq_w_rel_tol(m3d.get_U(0, j, k), mxz.get_U(j, k));
      if (!res) return res;
    }
  }

  return testing::AssertionSuccess();
}

template <class olim>
testing::AssertionResult
solution_is_exact_in_factored_square(
  int n, double tol = eps<double>, std::enable_if_t<olim::ndim == 2> * = 0)
{
  double h = 2./(n - 1);
  int i0 = n/2, j0 = n/2;
  typename olim::fac_src_t src {{(double) i0, (double) j0}, 1.0};
  olim o {{n, n}, h, s0<2>, {1., 1.}};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      o.set_fac_src({i, j}, &src);
    }
  }
  o.add_boundary_node({i0, j0});
  o.run();
  for (int i = 0; i < n; ++i) {
    double y = h*i - 1;
    for (int j = 0; j < n; ++j) {
      double x = h*j - 1;
      double u = std::hypot(x, y);
      double U = o.get_U({i, j});
      if (fabs(u - U) > tol*fabs(u) + tol) {
        return testing::AssertionFailure()
          << "|" << u << " - " << U << "| > " << tol << "*" << fabs(u)
          << " + " << tol;
      }
    }
  }
  return testing::AssertionSuccess();
}

template <class olim>
testing::AssertionResult
solution_is_exact_in_factored_square(
  int n, double tol = eps<double>, std::enable_if_t<olim::ndim == 3> * = 0)
{
  double h = 2./(n - 1);
  int i0 = n/2, j0 = n/2, k0 = n/2;

  typename olim::fac_src_t src {{(double) i0, (double) j0, (double) k0}, 1.0};
  olim o {{n, n, n}, h, s0<3>, {1., 1., 1.}};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        o.set_fac_src({i, j, k}, &src);
      }
    }
  }
  o.add_boundary_node({i0, j0, k0});
  o.run();
  for (int i = 0; i < n; ++i) {
    double x = h*i - 1;
    for (int j = 0; j < n; ++j) {
      double y = h*j - 1;
      for (int k = 0; k < n; ++k) {
        double z = h*k - 1;
        double u = std::sqrt(x*x + y*y + z*z);
        double U = o.get_U({i, j, k});
        if (fabs(u - U) > tol*fabs(u) + tol) {
          return testing::AssertionFailure()
            << "|" << u << " - " << U << "| > " << tol << "*" << fabs(u)
            << " + " << tol;
        }
      }
    }
  }
  return testing::AssertionSuccess();
}

#endif // __OLIM_TEST_COMMON_HPP__
