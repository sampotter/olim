#include <gtest/gtest.h>

#include "olim.hpp"

#include "olim_wrapper.h"

TEST (olim_wrapper, olim4_rhr_works) {
  int n = 11;
  double h = 1./(n - 1);
  int dims[2] = {6, 11};
  int nelts = dims[0]*dims[1];

  olim_wrapper_params p {OLIM4, RHR, h, 2, dims};
  olim_wrapper * w = nullptr;
  olim_wrapper_init(&w, &p);

  olim4_rhr olim {{dims}, h, no_slow_t {}};

  double * U = olim.get_U_ptr(), * U_wrap = nullptr;
  olim_wrapper_get_U_ptr(w, &U_wrap);
  ASSERT_NE(U_wrap, nullptr);
  for (int i = 0; i < nelts; ++i) {
    ASSERT_EQ(U[i], U_wrap[i]);
    ASSERT_EQ(U[i], inf<double>);
  }

  char * state = olim.get_state_ptr(), * state_wrap = nullptr;
  olim_wrapper_get_state_ptr(w, &state_wrap);
  ASSERT_NE(state_wrap, nullptr);
  for (int i = 0; i < nelts; ++i) {
    ASSERT_EQ(state[i], state_wrap[i]);
    ASSERT_EQ(state[i], static_cast<char>(state::far));
  }
}
