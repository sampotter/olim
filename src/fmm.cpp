#include "fmm.h"

#include <algorithm>

#include "basic_marcher.hpp"
#include "basic_marcher_3d.hpp"
#include "olim.hpp"
#include "olim3d.hpp"

void fmm(double * out, bool * in, int M, int N, double h, double * s_cache,
         marcher_type type) {
  std::unique_ptr<marcher<node>> m;
  if (s_cache == nullptr) {
    if (type == BASIC) {
      m = std::make_unique<basic_marcher>(M, N, h);
    } else if (type == OLIM4_MP0) {
      m = std::make_unique<olim4_mp0>(M, N, h);
    } else if (type == OLIM4_RHR) {
      m = std::make_unique<olim4_rhr>(M, N, h);
    } else if (type == OLIM8_MP0) {
      m = std::make_unique<olim8_mp0>(M, N, h);
    } else if (type == OLIM8_MP1) {
      m = std::make_unique<olim8_mp1>(M, N, h);
    } else if (type == OLIM8_RHR) {
      m = std::make_unique<olim8_rhr>(M, N, h);
    }
  } else {
    if (type == BASIC) {
      m = std::make_unique<basic_marcher>(M, N, h, s_cache);
    } else if (type == OLIM4_MP0) {
      m = std::make_unique<olim4_mp0>(M, N, h, s_cache);
    } else if (type == OLIM4_RHR) {
      m = std::make_unique<olim4_rhr>(M, N, h, s_cache);
    } else if (type == OLIM8_MP0) {
      m = std::make_unique<olim8_mp0>(M, N, h, s_cache);
    } else if (type == OLIM8_MP1) {
      m = std::make_unique<olim8_mp1>(M, N, h, s_cache);
    } else if (type == OLIM8_RHR) {
      m = std::make_unique<olim8_rhr>(M, N, h, s_cache);
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

void fmm3d(double * out, bool * in, int * dims, double h, double * s_cache,
           marcher_type type) {
  int M1 = dims[0], M2 = dims[1], M3 = dims[2];

  std::unique_ptr<marcher_3d<node_3d>> m;
  if (s_cache == nullptr) {
    if (type == BASIC) {
      m = std::make_unique<basic_marcher_3d>(M1, M2, M3, h);
    } else if (type == OLIM6_MP0) {
      m = std::make_unique<olim6_mp0>(M1, M2, M3, h);
    } else if (type == OLIM6_MP1) {
      m = std::make_unique<olim6_mp1>(M1, M2, M3, h);
    } else if (type == OLIM6_RHR) {
      m = std::make_unique<olim6_rhr>(M1, M2, M3, h);
    } else if (type == OLIM18_MP0) {
      m = std::make_unique<olim18_mp0>(M1, M2, M3, h);
    } else if (type == OLIM18_MP1) {
      m = std::make_unique<olim18_mp1>(M1, M2, M3, h);
    } else if (type == OLIM18_RHR) {
      m = std::make_unique<olim18_rhr>(M1, M2, M3, h);
    } else if (type == OLIM26_MP0) {
      m = std::make_unique<olim26_mp0>(M1, M2, M3, h);
    } else if (type == OLIM26_MP1) {
      m = std::make_unique<olim26_mp1>(M1, M2, M3, h);
    } else if (type == OLIM26_RHR) {
      m = std::make_unique<olim26_rhr>(M1, M2, M3, h);
    }
  } else {
    if (type == BASIC) {
      m = std::make_unique<basic_marcher_3d>(M1, M2, M3, h, s_cache);
    } else if (type == OLIM6_MP0) {
      m = std::make_unique<olim6_mp0>(M1, M2, M3, h, s_cache);
    } else if (type == OLIM6_MP1) {
      m = std::make_unique<olim6_mp1>(M1, M2, M3, h, s_cache);
    } else if (type == OLIM6_RHR) {
      m = std::make_unique<olim6_rhr>(M1, M2, M3, h, s_cache);
    } else if (type == OLIM18_MP0) {
      m = std::make_unique<olim18_mp0>(M1, M2, M3, h, s_cache);
    } else if (type == OLIM18_MP1) {
      m = std::make_unique<olim18_mp1>(M1, M2, M3, h, s_cache);
    } else if (type == OLIM18_RHR) {
      m = std::make_unique<olim18_rhr>(M1, M2, M3, h, s_cache);
    } else if (type == OLIM26_MP0) {
      m = std::make_unique<olim26_mp0>(M1, M2, M3, h, s_cache);
    } else if (type == OLIM26_MP1) {
      m = std::make_unique<olim26_mp1>(M1, M2, M3, h, s_cache);
    } else if (type == OLIM26_RHR) {
      m = std::make_unique<olim26_rhr>(M1, M2, M3, h, s_cache);
    }
  }

  for (int i = 0, l = 0; i < M1; ++i) {
    for (int j = 0; j < M2; ++j) {
      for (int k = 0; k < M3; ++k) {
        if (in[l++]) {
          m->add_boundary_node(i, j, k);
        }
      }
    }
  }

  m->run();

  for (int i = 0, l = 0; i < M1; ++i) {
    for (int j = 0; j < M2; ++j) {
      for (int k = 0; k < M3; ++k) {
        out[l++] = m->get_value(i, j, k);
      }
    }
  }
}
