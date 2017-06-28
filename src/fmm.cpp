#include "mex.h"

#include <cassert>
#include <string>

#include "fmm_mex.hpp"

static mxArray * getDefaultPvaluesMatrix(int M, int N) {
  mxArray * pvalues = mxCreateDoubleMatrix(M + 2, N + 2, mxREAL);
  double * pv = mxGetPr(pvalues);
  int nelts = (M + 2)*(N + 2);
  for (int i = 0; i < nelts; ++i) {
    pv[i] = 1;
  }
  return pvalues;
}

static double parseH(mxArray const * arg) {
  if (!mxIsClass(arg, "double") || mxGetM(arg) > 1 || mxGetN(arg) > 1) {
    mexErrMsgTxt("Argument containing h must be a double scalar.");
  }
  return mxGetScalar(arg);
}

static void parseSpeedFunc(mxArray * arg, mxArray *& pvalues, double h, int M,
                           int N) {
  if (!mxIsClass(arg, "function_handle")) {
    mexErrMsgTxt("Speed function argument must be a function.");
  }

  mxArray * X = mxCreateDoubleMatrix(M + 2, N + 2, mxREAL);
  mxArray * Y = mxCreateDoubleMatrix(M + 2, N + 2, mxREAL);

  double * Xpr = mxGetPr(X);
  double * Ypr = mxGetPr(Y);

  int k = 0;
  for (int i = -1; i <= M; ++i) {
    double y = h*i;
    for (int j = -1; j <= N; ++j) {
      Xpr[k] = h*j;
      Ypr[k++] = y;
    }
  }

  mxArray * plhs[] = {nullptr};
  mxArray * prhs[] = {arg, X, Y};
  mxArray * except = mexCallMATLABWithTrap(1, plhs, 3, prhs, "feval");
  assert(except == nullptr);
  mxDestroyArray(pvalues);
  pvalues = plhs[0];

  mxDestroyArray(X);
  mxDestroyArray(Y);
}

static marcher_type parseMarcherType(mxArray const * arg) {
  if (!mxIsClass(arg, "char") || mxGetM(arg) != 1) {
    mexErrMsgTxt("Argument containing method must be a string.");
  }
  std::string str {mxArrayToString(arg)};
  marcher_type type;
  if (str == "basic") {
    type = marcher_type::basic;
  } else if (str == "olim8pt") {
    type = marcher_type::olim8pt;
  } else {
    mexErrMsgTxt(("Invalid marcher type: " + str).c_str());
  }
  return type;
}

static void
parseKeywordArguments(int nlhs, mxArray * plhs[],
                      int nrhs, mxArray const * prhs[],
                      marcher_type & type, double & h, mxArray *& pvalues,
                      int M, int N) {
  if (nrhs % 2 == 0) {
    mexErrMsgTxt("Keyword arguments must be passed in pairs.");
  }
  // Need to parse h first to use later
  for (int i = 1; i < nrhs; i += 2) {
    if (!mxIsClass(prhs[i], "char") || mxGetM(prhs[i]) != 1) {
      mexErrMsgTxt("Keywords must be strings.");
    }
    if (std::string {mxArrayToString(prhs[i])} == "h") {
      h = parseH(prhs[i + 1]);
    }
  }
  for (int i = 1; i < nrhs; i += 2) {
    if (!mxIsClass(prhs[i], "char") || mxGetM(prhs[i]) != 1) {
      mexErrMsgTxt("Keywords must be strings.");
    }
    std::string keyword {mxArrayToString(prhs[i])};
    if (keyword == "h") {
      continue;
    }
    if (keyword == "Speed") {
      parseSpeedFunc(const_cast<mxArray *>(prhs[i + 1]), pvalues, h, M, N);
    } else if (keyword == "Method") {
      type = parseMarcherType(prhs[i + 1]);
    } else {
      mexErrMsgTxt(("Invalid keyword: " + keyword).c_str());
    }
  }
}

static void
parseRegularArguments(int nlhs, mxArray * plhs[],
                      int nrhs, mxArray const * prhs[],
                      marcher_type & type, double & h, mxArray *& pvalues,
                      int M, int N) {
  if (nrhs >= 2) {
    h = parseH(prhs[1]);
  }
  if (nrhs >= 3) {
    parseSpeedFunc(const_cast<mxArray *>(prhs[2]), pvalues, h, M, N);
  }
  if (nrhs >= 4) {
    type = parseMarcherType(prhs[3]);
  }
}

/**
 * nlhs -- number of outputs (expected: 1)
 * plhs -- outputs: a double matrix
 * nrhs -- number of inputs (expected: >=1)
 * prhs -- inputs: a logical matrix (the boundary points), a double
 *         (the uniform grid spacing), a function handle or matrix (the
 *         speed function), a string (which method to use)
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[]) {
  /**
   * Ensure the correct number of outputs.
   */
  if (nlhs > 1) {
    mexErrMsgTxt("Too many outputs.");
  }

  /**
   * Get the height (M) and width (N) of the input boundary matrix.
   */
  int M = mxGetM(prhs[0]);
  int N = mxGetN(prhs[0]);

  /**
   * Set up arguments from passed parameters.
   */
  marcher_type type = marcher_type::basic;
  double h = 1;
  mxArray * pvalues = getDefaultPvaluesMatrix(M, N);

  if (nrhs >= 2 && mxIsClass(prhs[1], "char") && mxGetM(prhs[1]) == 1) {
    parseKeywordArguments(nlhs, plhs, nrhs, prhs, type, h, pvalues, M, N);
  } else {
    parseRegularArguments(nlhs, plhs, nrhs, prhs, type, h, pvalues, M, N);
  }

  /**
   * Set up input arguments on our end.
   */
  mxLogical * logicals = mxGetLogicals(prhs[0]);
  bool * in = (bool *) calloc(M*N, sizeof(bool));
  for (int i = 0; i < M*N; ++i) {
    in[i] = logicals[i];
  }

  /**
   * Initialize output and get pointer.
   */
  plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
  double * out = mxGetPr(plhs[0]);

  fmm_mex(out, in, M, N, h, mxGetPr(pvalues), type);

  /**
   * Clean up.
   */
  mxDestroyArray(pvalues);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
