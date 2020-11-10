#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <eikonal/olim3d.hpp>

constexpr ordering ord = ordering::COLUMN_MAJOR;

using ivec = vec<int, 3>;

bool is_space(char c) {
  return c == ' ' || c == '\t';
}

bool is_eol(char c) {
  return c == '\n' || c == '\0';
}

bool is_sep(char c) {
  return is_space(c) || is_eol(c);
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("usage: %s <path>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  FILE *fp = nullptr;
  fp = fopen(argv[1], "rb");

  ivec dim;

  char buf[64];

  char c;
  for (int m = 0; m < 3; ++m) {
    int i = 0;
    while (!is_sep(c = fgetc(fp))) {
      buf[i++] = c;
    }
    buf[i] = '\0';
    dim[m] = atoi(buf);
  }
  while (!is_eol(c)) {
    c = fgetc(fp);
  }

  printf("shape: %d, %d, %d\n", dim[0], dim[1], dim[2]);

  size_t nelts = dim[0]*dim[1]*dim[2];
  float *onsets = (float *)malloc(sizeof(float)*nelts);
  size_t nread = fread(onsets, sizeof(float), nelts, fp);
  fclose(fp);
  if (nelts != nread) {
    printf("only read %lu items (expected %lu)\n", nread, nelts);
    free(onsets);
    exit(EXIT_FAILURE);
  }

  double *s = (double *)malloc(sizeof(double)*nelts);
  for (size_t l = 0; l < nelts; ++l) {
    s[l] = 1.;
  }

  double h = 0.015937;
  eikonal::olim3d_rhr<>::ivec dims(dim[0], dim[1], dim[2]);
  eikonal::olim3d_rhr<> olim(dims, h, s);

  float min_onset = INFINITY;
  ivec ind, ind_src;
  for (size_t l = 0; l < nelts; ++l) {
    if (std::isnan(onsets[l])) {
      ind = to_vector_index<ord>(l, dim);
      olim.add_bd(ind);
    } else if (onsets[l] < min_onset) {
      ind_src = to_vector_index<ord>(l, dim);
      min_onset = onsets[l];
    }
  }

  printf("T(%d, %d, %d) = %f\n", ind_src[0], ind_src[1], ind_src[2], min_onset);

  size_t l = to_linear_index<ord>(ind, dim);
  printf("T[%lu] = %f\n", l, onsets[l]);

  olim.add_src(ind_src);
  olim.solve();

  free(onsets);

  double *T = olim.get_U_ptr();

  size_t nelts_padded = (dims[0] + 2)*(dims[1] + 2)*(dims[2] + 2);

  fp = fopen("T.bin", "wb");
  fwrite(T, sizeof(double), nelts_padded, fp);
  fclose(fp);

  // Quick groundtruth test

  eikonal::olim3d_rhr<> olim_(dims, h, s);
  olim_.add_src(ind_src);
  olim_.solve();
  T = olim_.get_U_ptr();

  fp = fopen("T_gt.bin", "wb");
  fwrite(T, sizeof(double), nelts_padded, fp);
  fclose(fp);

  free(s);
}
