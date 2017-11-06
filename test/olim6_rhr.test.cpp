#include "test.hpp"

#include "basic_marcher_3d.hpp"
#include "olim.test.common.hpp"
#include "olim4.hpp"
#include "olim6.hpp"

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

void agrees_with_basic_marcher_3d() {
  int n = 11;
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
  int n = 11;
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
  using olim = olim4_rhr;
  using olim3d = olim6_rhr;

  quadrants_are_correct<olim3d>(1 + sqrt(2)/2);
  octants_are_correct<olim3d>(
    1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3);
  planes_are_correct<olim, olim3d>();
  planes_are_correct_for_nontrivial_speed_function();
  result_is_symmetric<olim3d>();
  result_is_symmetric_for_nontrivial_speed_function();
  two_by_two_by_three_cells_are_correct<olim3d>();
  agrees_with_basic_marcher_3d();
  agrees_with_basic_marcher_3d_for_nontrivial_speed_function();
  plane_boundaries_are_correct();
}
