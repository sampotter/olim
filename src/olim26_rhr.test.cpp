#include "test.hpp"

#include <cmath>

#include "olim8_rhr.hpp"
#include "olim26_rhr.hpp"


void quadrants_are_correct() {
  int n = 2;
  double h = 1;

  /*
   * Tests for quadrants in the x-y plane:
   */
  {
    olim26_rhr_arma m {n, n, 1, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), sqrt(2));
  }
  {
    olim26_rhr_arma m {n, n, 1, h, default_speed_func_3d, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim26_rhr_arma m {n, n, 1, h, default_speed_func_3d, 0, 1, 0};
    m.add_boundary_node(1, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim26_rhr_arma m {n, n, 1, h, default_speed_func_3d, 1, 1, 0};
    m.add_boundary_node(1, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 0.0);
  }

  /**
   * Tests for quadrants in the x-z plane:
   */
  {
    olim26_rhr_arma m {n, 1, n, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), sqrt(2));
  }
  {
    olim26_rhr_arma m {n, 1, n, h, default_speed_func_3d, 0, 1, 0};
    m.add_boundary_node(1, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0);
  }
  {
    olim26_rhr_arma m {n, 1, n, h, default_speed_func_3d, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0);
  }
  {
    olim26_rhr_arma m {n, 1, n, h, default_speed_func_3d, 0, 1, 1};
    m.add_boundary_node(1, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 0.0);
  }

  /**
   * Tests for quadrants in the y-z plane:
   */
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), sqrt(2));
  }
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 1, 0, 1};
    m.add_boundary_node(0, 1, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 0.0);
  }
}

void two_by_two_by_two_octants_are_correct() {
  int n = 2;
  double h = 1;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        double x0 = j, y0 = i, z0 = k;
        olim26_rhr_arma m {n, n, n, h, default_speed_func_3d, x0, y0, z0};
        m.add_boundary_node(i, j, k);
        m.run();
        IS_APPROX_EQUAL(m.get_value(i, j, k), 0.0);
        IS_APPROX_EQUAL(m.get_value((i + 1) % 2, j, k), 1.0);
        IS_APPROX_EQUAL(m.get_value(i, (j + 1) % 2, k), 1.0);
        IS_APPROX_EQUAL(m.get_value(i, j, (k + 1) % 2), 1.0);
        IS_APPROX_EQUAL(m.get_value((i + 1) % 2, (j + 1) % 2, k), sqrt(2));
        IS_APPROX_EQUAL(m.get_value((i + 1) % 2, j, (k + 1) % 2), sqrt(2));
        IS_APPROX_EQUAL(m.get_value(i, (j + 1) % 2, (k + 1) % 2), sqrt(2));
        IS_APPROX_EQUAL(
          m.get_value((i + 1) % 2, (j + 1) % 2, (k + 1) % 2), sqrt(3));
      }
    }
  }
}

void planes_are_correct() {
  int n = 3;
  double h = 1.0/(n/2);
  
  olim8_rhr m8 {n, n, h, default_speed_func, 1, 1};
  m8.add_boundary_node(n/2, n/2);
  m8.run();
  
  olim26_rhr_arma m26 {n, n, n, h, default_speed_func_3d, 1, 1, 1};
  m26.add_boundary_node(n/2, n/2, n/2);
  m26.run();

  /**
   * Tests for quadrants in the y-z plane:
   */
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), sqrt(2));
  }
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 1, 0, 1};
    m.add_boundary_node(0, 1, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 0.0);
  }
}

void result_is_symmetric() {
  int n = 5;
  double x0 = (n - 1.0)/2.0, y0 = x0, z0 = x0;
  olim26_rhr_arma m {n, n, n, 1.0, default_speed_func_3d, x0, y0, z0};
  m.add_boundary_node(n/2, n/2, n/2);
  m.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0, k_ = n - 1; k < n; ++k, --k_) {
        IS_APPROX_EQUAL(m.get_value(i, j, k), m.get_value(i, j, k_));
      }
    }
  }
}

void two_by_two_by_three_cells_are_correct() {
  int dims[3][3] = {{3, 2, 2}, {2, 3, 2}, {2, 2, 3}};

  olim26_rhr_arma m_gt {3, 2, 2, 1, default_speed_func_3d, 0, 0, 0};
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
      olim26_rhr_arma m {imax, jmax, kmax, 1, default_speed_func_3d, x0, y0, z0};
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

void plane_boundaries_are_correct() {
  int n = 2;
  double h = 1;
  olim26_rhr_arma m {n, n, n, h, default_speed_func_3d, 0, 0, 0};
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
  two_by_two_by_two_octants_are_correct();
  planes_are_correct();
  result_is_symmetric();
  two_by_two_by_three_cells_are_correct();
  plane_boundaries_are_correct();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
