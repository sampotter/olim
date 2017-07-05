#include "fmm_mex.hpp"

#include <memory>

#include "basic_marcher.hpp"
#include "olim8_mp0.hpp"
#include "olim8_mp1.hpp"
#include "olim8_rhr.hpp"

void fmm_mex(double * out, bool * in, int M, int N, double h,
			 double * S, marcher_type type) {
  std::unique_ptr<fast_marcher> m;
  if (type == marcher_type::basic) {
    m = std::make_unique<basic_marcher>(M, N, h, S);
  } else if (type == marcher_type::olim8_mp0) {
    m = std::make_unique<olim8_mp0>(M, N, h, S);
  } else if (type == marcher_type::olim8_mp1) {
    m = std::make_unique<olim8_mp1>(M, N, h, S);
  } else if (type == marcher_type::olim8_rhr) {
    m = std::make_unique<olim8_rhr>(M, N, h, S);
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
