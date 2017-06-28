#include "fmm_mex.hpp"

#include <memory>

#include "basic_marcher.hpp"
#include "olim_8pt_rhr.hpp"

void fmm_mex(double * out, bool * in, size_t M, size_t N, double h,
			 double * S, marcher_type type) {
  std::unique_ptr<fast_marcher> m;
  if (type == marcher_type::basic) {
    m = std::make_unique<basic_marcher>(M, N, h, S);
  } else if (type == marcher_type::olim8pt) {
    m = std::make_unique<olim_8pt_rhr>(M, N, h, S);
  }
	
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (in[N*i + j]) {
        m->add_boundary_node(i, j);
      }
    }
  }

  m->run();

  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      out[M*i + j] = m->get_value(i, j);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
