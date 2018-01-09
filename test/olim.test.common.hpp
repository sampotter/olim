#ifndef __OLIM_TEST_COMMON_HPP__
#define __OLIM_TEST_COMMON_HPP__

#include <cassert>
#include <vector>

#include "speed_funcs.hpp"
#include "test.hpp"
#include "typedefs.h"

template <class olim>
void trivial_case_works() {
  olim m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, 0), 0.0);
}

template <class olim>
void adjacent_update_works() {
  olim m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(1, 0), 0.5);
}

template <class olim>
void correct_corners_in_limit(int n, double tol) {
  double h = 2./(n - 1);
  olim m {n, n, h, (speed_func) default_speed_func, 1., 1.};
  m.add_boundary_node(n/2, n/2);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, n - 1), sqrt(2), tol);
  IS_APPROX_EQUAL(m.get_value(0, n - 1), sqrt(2), tol);
  IS_APPROX_EQUAL(m.get_value(n - 1, 0), sqrt(2), tol);
  IS_APPROX_EQUAL(m.get_value(n - 1, n - 1), sqrt(2), tol);
}

template <class olim3d>
void quadrants_are_correct(double diag_value) {
  int n = 2;
  double h = 1;

  /*
   * Tests for quadrants in the x-y plane:
   */
  {
    olim3d m {n, n, 1, h, (speed_func_3d) default_speed_func, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), diag_value);
  }
  {
    olim3d m {n, n, 1, h, (speed_func_3d) default_speed_func, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), diag_value);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim3d m {n, n, 1, h, (speed_func_3d) default_speed_func, 0, 1, 0};
    m.add_boundary_node(1, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), diag_value);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim3d m {n, n, 1, h, (speed_func_3d) default_speed_func, 1, 1, 0};
    m.add_boundary_node(1, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), diag_value);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 0.0);
  }

  /**
   * Tests for quadrants in the x-z plane:
   */
  {
    olim3d m {n, 1, n, h, (speed_func_3d) default_speed_func, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), diag_value);
  }
  {
    olim3d m {n, 1, n, h, (speed_func_3d) default_speed_func, 0, 1, 0};
    m.add_boundary_node(1, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), diag_value);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0);
  }
  {
    olim3d m {n, 1, n, h, (speed_func_3d) default_speed_func, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), diag_value);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0);
  }
  {
    olim3d m {n, 1, n, h, (speed_func_3d) default_speed_func, 0, 1, 1};
    m.add_boundary_node(1, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), diag_value);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 0.0);
  }

  /**
   * Tests for quadrants in the y-z plane:
   */
  {
    olim3d m {1, n, n, h, (speed_func_3d) default_speed_func, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), diag_value);
  }
  {
    olim3d m {1, n, n, h, (speed_func_3d) default_speed_func, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), diag_value);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim3d m {1, n, n, h, (speed_func_3d) default_speed_func, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), diag_value);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim3d m {1, n, n, h, (speed_func_3d) default_speed_func, 1, 0, 1};
    m.add_boundary_node(0, 1, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), diag_value);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 0.0);
  }
}

template <class olim3d>
void octants_are_correct(double diag2val, double diag3val) {
  int n = 2;
  double h = 1;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        double x0 = j, y0 = i, z0 = k;
        olim3d m {n, n, n, h, (speed_func_3d) default_speed_func, x0, y0, z0};
        m.add_boundary_node(i, j, k);
        m.run();
        IS_APPROX_EQUAL(m.get_value(i, j, k), 0.0);
        IS_APPROX_EQUAL(m.get_value((i + 1) % 2, j, k), 1.0);
        IS_APPROX_EQUAL(m.get_value(i, (j + 1) % 2, k), 1.0);
        IS_APPROX_EQUAL(m.get_value(i, j, (k + 1) % 2), 1.0);
        IS_APPROX_EQUAL(m.get_value((i + 1) % 2, (j + 1) % 2, k), diag2val);
        IS_APPROX_EQUAL(m.get_value((i + 1) % 2, j, (k + 1) % 2), diag2val);
        IS_APPROX_EQUAL(m.get_value(i, (j + 1) % 2, (k + 1) % 2), diag2val);
        IS_APPROX_EQUAL(m.get_value((i + 1) % 2, (j + 1) % 2, (k + 1) % 2),
                        diag3val);
      }
    }
  }
}

template <class olim, class olim3d>
void planes_are_correct(
  speed_func s = (speed_func) default_speed_func,
  speed_func_3d s3d = (speed_func_3d) default_speed_func,
  int n = 5)
{
  assert(n % 2 == 1);

  double h = 2.0/(n - 1);
  
  olim m {n, n, h, s, 1, 1};
  m.add_boundary_node(n/2, n/2);
  m.run();
  
  olim3d m3d {n, n, n, h, s3d, 1, 1, 1};
  m3d.add_boundary_node(n/2, n/2, n/2);
  m3d.run();

  // Check that planes are correct:
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double U = m.get_value(i, j);
      IS_APPROX_EQUAL(U, m3d.get_value(i, j, n/2));
      IS_APPROX_EQUAL(U, m3d.get_value(i, n/2, j));
      IS_APPROX_EQUAL(U, m3d.get_value(n/2, i, j));
    }
  }
}

template <class olim>
void result_is_symmetric(speed_func s = default_speed_func, int n = 11) {
  double x0 = (n - 1)/2.0;
  olim m {n, n, 1.0, s, x0, x0};
  m.add_boundary_node(n/2, n/2);
  m.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0, j_ = n - 1; j < n; ++j, --j_) {
      IS_APPROX_EQUAL(m.get_value(i, j), m.get_value(i, j_));
    }
  }

  for (int i = 0, i_ = n - 1; i < n; ++i, --i_) {
    for (int j = 0; j < n; ++j) {
      IS_APPROX_EQUAL(m.get_value(i, j), m.get_value(i_, j));
    }
  }
}

template <class olim3d>
void result_is_symmetric(speed_func_3d s = default_speed_func, int n = 5) {
  double x0 = (n - 1)/2.0;
  olim3d m {n, n, n, 1.0, s, x0, x0, x0};
  m.add_boundary_node(n/2, n/2, n/2);
  m.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0, k_ = n - 1; k < n; ++k, --k_) {
        IS_APPROX_EQUAL(m.get_value(i, j, k), m.get_value(i, j, k_));
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      for (int j = 0, j_ = n - 1; j < n; ++j, --j_) {
        IS_APPROX_EQUAL(m.get_value(i, j, k), m.get_value(i, j_, k));
      }
    }
  }

  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) {
      for (int i = 0, i_ = n - 1; i < n; ++i, --i_) {
        IS_APPROX_EQUAL(m.get_value(i, j, k), m.get_value(i_, j, k));
      }
    }
  }
}

template <class olim3d>
void two_by_two_by_three_cells_are_correct() {
  int dims[3][3] = {{3, 2, 2}, {2, 3, 2}, {2, 2, 3}};

  olim3d m_gt {3, 2, 2, 1, (speed_func_3d) default_speed_func, 0, 0, 0};
  m_gt.add_boundary_node(0, 0, 0);
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
      olim3d m {imax, jmax, kmax, 1, (speed_func_3d) default_speed_func, x0, y0, z0};
      m.add_boundary_node(i0, j0, k0);
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
            IS_APPROX_EQUAL(m.get_value(i, j, k), m_gt.get_value(a_, b_, c_));
          }
        }
      }
    }
  }
}

template <class olim3d>
void plane_boundaries_are_correct() {
  int n = 2;
  double h = 1;
  olim3d m {n, n, n, h, (speed_func_3d) default_speed_func, 0, 0, 0};

  typename olim3d::node_type nodes[4];
  for (int i = 0, k = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      nodes[k++] = typename olim3d::node_type {i, j, 0};
    }
  }
  m.add_boundary_nodes(nodes, 4);

  m.run();
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      IS_APPROX_EQUAL(m.get_value(i, j, 1), 1.0);
    }
  }
}

template <class olim1, class olim2>
void olims_agree(speed_func s = default_speed_func, int n = 11) {
  double h = 2.0/(n - 1);
  int i0 = n/2;

  olim1 m1 {n, n, h, s, 1, 1};
  m1.add_boundary_node(i0, i0);
  m1.run();

  olim2 m2 {n, n, h, s, 1, 1};
  m2.add_boundary_node(i0, i0);
  m2.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      IS_APPROX_EQUAL(m1.get_value(i, j), m2.get_value(i, j));
    }
  }
}

template <class olim3d, class other_olim3d>
void agrees_with_other_olim3d(int n = 11) {
  double h = 2.0/(n - 1);
  int i0 = n/2;

  olim3d m1 {n, n, n, h, (speed_func_3d) default_speed_func, 1, 1, 1};
  m1.add_boundary_node(i0, i0, i0);
  m1.run();

  other_olim3d m2 {n, n, n, h, (speed_func_3d) default_speed_func, 1, 1, 1};
  m2.add_boundary_node(i0, i0, i0);
  m2.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        IS_APPROX_EQUAL(m1.get_value(i, j, k), m2.get_value(i, j, k));
      }
    }
  }
}

#endif // __OLIM_TEST_COMMON_HPP__
