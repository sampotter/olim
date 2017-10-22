#include "test.hpp"

#include "basic_marcher_3d.hpp"
#include "olim4_rhr.hpp"
#include "olim6.hpp"
#include "speed_funcs.hpp"

void quadrants_are_correct() {
  int n = 2;
  double h = 1;

  /*
   * Tests for quadrants in the x-y plane:
   */
  {
    olim6_rhr m {n, n, 1, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0 + sqrt(2)/2.0);
  }
  {
    olim6_rhr m {n, n, 1, h, default_speed_func_3d, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim6_rhr m {n, n, 1, h, default_speed_func_3d, 0, 1, 0};
    m.add_boundary_node(1, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim6_rhr m {n, n, 1, h, default_speed_func_3d, 1, 1, 0};
    m.add_boundary_node(1, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 0.0);
  }

  /**
   * Tests for quadrants in the x-z plane:
   */
  {
    olim6_rhr m {n, 1, n, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0 + sqrt(2)/2.0);
  }
  {
    olim6_rhr m {n, 1, n, h, default_speed_func_3d, 0, 1, 0};
    m.add_boundary_node(1, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0);
  }
  {
    olim6_rhr m {n, 1, n, h, default_speed_func_3d, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0);
  }
  {
    olim6_rhr m {n, 1, n, h, default_speed_func_3d, 0, 1, 1};
    m.add_boundary_node(1, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 0.0);
  }

  /**
   * Tests for quadrants in the y-z plane:
   */
  {
    olim6_rhr m {1, n, n, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0 + sqrt(2)/2.0);
  }
  {
    olim6_rhr m {1, n, n, h, default_speed_func_3d, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim6_rhr m {1, n, n, h, default_speed_func_3d, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim6_rhr m {1, n, n, h, default_speed_func_3d, 1, 0, 1};
    m.add_boundary_node(0, 1, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 0.0);
  }
}

void octants_are_correct() {
  int n = 2;
  double h = 1;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        double x0 = j, y0 = i, z0 = k;
        olim6_rhr m {n, n, n, h, default_speed_func_3d, x0, y0, z0};
        m.add_boundary_node(i, j, k);
        m.run();
        IS_APPROX_EQUAL(m.get_value(i, j, k), 0.0);
        IS_APPROX_EQUAL(m.get_value((i + 1) % 2, j, k), 1.0);
        IS_APPROX_EQUAL(m.get_value(i, (j + 1) % 2, k), 1.0);
        IS_APPROX_EQUAL(m.get_value(i, j, (k + 1) % 2), 1.0);
        IS_APPROX_EQUAL(
          m.get_value((i + 1) % 2, (j + 1) % 2, k), 1.0 + sqrt(2)/2.0);
        IS_APPROX_EQUAL(
          m.get_value((i + 1) % 2, j, (k + 1) % 2), 1.0 + sqrt(2)/2.0);
        IS_APPROX_EQUAL(
          m.get_value(i, (j + 1) % 2, (k + 1) % 2), 1.0 + sqrt(2)/2.0);
        IS_APPROX_EQUAL(
          m.get_value((i + 1) % 2, (j + 1) % 2, (k + 1) % 2),
          1.0 + sqrt(2)/2.0 + sqrt(3)/3.0);
      }
    }
  }
}

void planes_are_correct() {
  int n = 11;
  double h = 1.0/(n/2);
  
  olim4_rhr m4 {n, n, h, default_speed_func, 1, 1};
  m4.add_boundary_node(n/2, n/2);
  m4.run();
  
  olim6_rhr m6 {n, n, n, h, default_speed_func_3d, 1, 1, 1};
  m6.add_boundary_node(n/2, n/2, n/2);
  m6.run();

  // Check that planes are correct:
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double m4value = m4.get_value(i, j);
      IS_APPROX_EQUAL(m4value, m6.get_value(i, j, n/2));
      IS_APPROX_EQUAL(m4value, m6.get_value(i, n/2, j));
      IS_APPROX_EQUAL(m4value, m6.get_value(n/2, i, j));
    }
  }
}

void planes_are_correct_for_nontrivial_speed_function() {
  int n = 11;
  double h = 1.0/(n/2);

  olim4_rhr m4 {n, n, h, (speed_func) s1, 1, 1};
  m4.add_boundary_node(n/2, n/2);
  m4.run();

  olim6_rhr m6 {n, n, n, h, (speed_func_3d) s1, 1, 1, 1};
  m6.add_boundary_node(n/2, n/2, n/2);
  m6.run();

  // Check that planes are correct:
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double m4value = m4.get_value(i, j);
      IS_APPROX_EQUAL(m4value, m6.get_value(i, j, n/2));
      IS_APPROX_EQUAL(m4value, m6.get_value(i, n/2, j));
      IS_APPROX_EQUAL(m4value, m6.get_value(n/2, i, j));
    }
  }
}

void result_is_symmetric() {
  int n = 5;
  olim6_rhr m {n, n, n, 1.0, default_speed_func_3d, 1.0, 1.0, 1.0};
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

void result_is_symmetric_for_nontrivial_speed_function() {
  int n = 5;
  olim6_rhr m {n, n, n, 0.5, (speed_func_3d) s1, 1.0, 1.0, 1.0};
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

void two_by_two_by_three_cells_are_correct() {
  int dims[3][3] = {{3, 2, 2}, {2, 3, 2}, {2, 2, 3}};

  basic_marcher_3d m_gt {3, 2, 2, 1, default_speed_func_3d, 0, 0, 0};
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
      olim6_rhr m {imax, jmax, kmax, 1, default_speed_func_3d, x0, y0, z0};
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

void agrees_with_basic_marcher_3d() {
  int n = 21;
  double h = 2.0/(n - 1);
  int i0 = (n - 1)/2, j0 = i0, k0 = i0;

  basic_marcher_3d m3d {n, n, n, h, default_speed_func_3d, 1, 1, 1};
  m3d.add_boundary_node(i0, j0, k0);
  m3d.run();

  olim6_rhr m6 {n, n, n, h, default_speed_func_3d, 1, 1, 1};
  m6.add_boundary_node(i0, j0, k0);
  m6.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        IS_APPROX_EQUAL(m3d.get_value(i, j, k), m6.get_value(i, j, k), 1e-13);
      }
    }
  }
}

void agrees_with_basic_marcher_3d_for_nontrivial_speed_function() {
  int n = 15;
  double h = 2.0/(n - 1);
  int i0 = (n - 1)/2, j0 = i0, k0 = i0;

  basic_marcher_3d m3d {n, n, n, h, (speed_func_3d) s1, 1, 1, 1};
  m3d.add_boundary_node(i0, j0, k0);
  m3d.run();

  olim6_rhr m6 {n, n, n, h, (speed_func_3d) s1, 1, 1, 1};
  m6.add_boundary_node(i0, j0, k0);
  m6.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        IS_APPROX_EQUAL(m3d.get_value(i, j, k), m6.get_value(i, j, k), 1e-13);
      }
    }
  }
}

void plane_boundaries_are_correct() {
  int n = 2;
  double h = 1;
  olim6_rhr m {n, n, n, h, default_speed_func_3d, 0, 0, 0};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      m.add_boundary_node(i, j, 0);
    }
  }
  m.run();
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      IS_APPROX_EQUAL(m.get_value(i, j, 1), 1.0);
    }
  }
}

int main() {
  quadrants_are_correct();
  octants_are_correct();
  planes_are_correct();
  planes_are_correct_for_nontrivial_speed_function();
  result_is_symmetric();
  result_is_symmetric_for_nontrivial_speed_function();
  two_by_two_by_three_cells_are_correct();
  agrees_with_basic_marcher_3d();
  agrees_with_basic_marcher_3d_for_nontrivial_speed_function();
  plane_boundaries_are_correct();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
