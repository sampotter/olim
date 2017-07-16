#include "fmm.h"

#include <memory>

#include "basic_marcher.hpp"
#include "olim4_mp0.hpp"
#include "olim4_rhr.hpp"
#include "olim4_rhr_lut.hpp"
#include "olim8_mp0.hpp"
#include "olim8_mp1.hpp"
#include "olim8_rhr.hpp"

void fmm(double * out, bool * in, int M, int N, double h, double * S,
         marcher_type type) {
  std::unique_ptr<marcher_2d> m;
  if (S == nullptr) {
    if (type == BASIC) {
      m = std::make_unique<basic_marcher>(M, N, h);
    } else if (type == OLIM4_MP0) {
      m = std::make_unique<olim4_mp0>(M, N, h);
    } else if (type == OLIM4_RHR) {
      m = std::make_unique<olim4_rhr>(M, N, h);
    } else if (type == OLIM4_RHR_LUT) {
      m = std::make_unique<olim4_rhr_lut>(M, N, h);
    } else if (type == OLIM8_MP0) {
      m = std::make_unique<olim8_mp0>(M, N, h);
    } else if (type == OLIM8_MP1) {
      m = std::make_unique<olim8_mp1>(M, N, h);
    } else if (type == OLIM8_RHR) {
      m = std::make_unique<olim8_rhr>(M, N, h);
    }
  } else {
    if (type == BASIC) {
      m = std::make_unique<basic_marcher>(M, N, h, S);
    } else if (type == OLIM4_MP0) {
      m = std::make_unique<olim4_mp0>(M, N, h, S);
    } else if (type == OLIM4_RHR) {
      m = std::make_unique<olim4_rhr>(M, N, h, S);
    } else if (type == OLIM4_RHR_LUT) {
      m = std::make_unique<olim4_rhr_lut>(M, N, h, S);
    } else if (type == OLIM8_MP0) {
      m = std::make_unique<olim8_mp0>(M, N, h, S);
    } else if (type == OLIM8_MP1) {
      m = std::make_unique<olim8_mp1>(M, N, h, S);
    } else if (type == OLIM8_RHR) {
      m = std::make_unique<olim8_rhr>(M, N, h, S);
    }
  }
	
  for (int i = 0, k = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (in[k++]) m->add_boundary_node(i, j);
    }
  }

  m->run();

  for (int j = 0, k = 0; j < N; ++j) {
    for (int i = 0; i < M; ++i) {
      out[k++] = m->get_value(i, j);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
