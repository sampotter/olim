#include "basic_marcher.hpp"
#include "basic_marcher_3d.hpp"
#include "olim.test.common.hpp"
#include "test.hpp"

void result_is_symmetric() {
  int n = 5;
  basic_marcher_3d m {n, n, n, 1.0, default_speed_func_3d, 1.0, 1.0, 1.0};
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
      basic_marcher_3d m {imax, jmax, kmax, 1, default_speed_func_3d, x0, y0, z0};
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
  basic_marcher_3d m {n, n, n, h, default_speed_func_3d, 0, 0, 0};

  node_3d nodes[4];
  for (int i = 0, k = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      nodes[k++] = node_3d {i, j, 0};
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

int main() {
  quadrants_are_correct<basic_marcher_3d>(1.0 + sqrt(2)/2);
  octants_are_correct<basic_marcher_3d>(
    1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3);
  planes_are_correct<basic_marcher, basic_marcher_3d>();
  result_is_symmetric();
  two_by_two_by_three_cells_are_correct();
  plane_boundaries_are_correct();
}
