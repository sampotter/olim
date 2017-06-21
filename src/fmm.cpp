#include "mex.h"

#include <string>

#include "fmm_mex.hpp"

typedef double (*speed_func)(double, double);

static double defaultSpeedFunc(double x, double y) {
  (void) x;
  (void) y;
  return 1;
}

static mxArray* userFuncPlhs[1] = {NULL};
static mxArray* userFuncPrhs[3] = {NULL, NULL, NULL};
static double* userFuncArgs = NULL;

static void setScalar(mxArray * arg, double value) {
  double * pr = mxGetPr(arg);
  *pr = value;
}

static double userSpeedFunc(double x, double y) {
  setScalar(userFuncPrhs[1], x);
  setScalar(userFuncPrhs[2], y);
  if (!mexCallMATLAB(1, userFuncPlhs, 3, userFuncPrhs, "feval")) {
    // TODO: error handling
  }
  return mxGetScalar(userFuncPlhs[0]);
}

static double userSlownessFunc(double x, double y) {
  setScalar(userFuncPrhs[1], x);
  setScalar(userFuncPrhs[2], y);
  if (!mexCallMATLAB(1, userFuncPlhs, 3, userFuncPrhs, "feval")) {
    // TODO: error handling
  }
  return 1.0/mxGetScalar(userFuncPlhs[0]);
}

static double parseH(mxArray const * arg) {
  if (!mxIsClass(arg, "double") || mxGetM(arg) > 1 || mxGetN(arg) > 1) {
    mexErrMsgTxt("Argument containing h must be a double scalar.");
  }
  return mxGetScalar(arg);
}

static void initUserFunc(mxArray const * arg) {
  userFuncPlhs[0] = mxCreateDoubleScalar(0);
  userFuncPrhs[0] = const_cast<mxArray *>(arg);
  userFuncPrhs[1] = mxCreateDoubleScalar(0);
  userFuncPrhs[2] = mxCreateDoubleScalar(0);
  userFuncArgs = static_cast<double *>(mxCalloc(2, sizeof(double)));
}

static void cleanupUserFunc() {
  mxDestroyArray(userFuncPlhs[0]);
  mxDestroyArray(userFuncPrhs[1]);
  mxDestroyArray(userFuncPrhs[2]);
  mxFree(userFuncArgs);
}

static speed_func parseSpeedFunc(mxArray const * arg) {
  if (!mxIsClass(arg, "function_handle")) {
    mexErrMsgTxt("Speed function argument is not a function handle.");
  }
  initUserFunc(arg);
  return userSpeedFunc;
}

static speed_func parseSlownessFunc(mxArray const * arg) {
  if (!mxIsClass(arg, "function_handle")) {
    mexErrMsgTxt("Slowness function argument is not a function handle.");
  }
  initUserFunc(arg);
  return userSlownessFunc;
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
                      marcher_type & type, double & h, speed_func & F) {
  if (nrhs % 2 == 0) {
    mexErrMsgTxt("Keyword arguments must be passed in pairs.");
  }
  for (int i = 1; i < nrhs; i += 2) {
    if (!mxIsClass(prhs[i], "char") || mxGetM(prhs[i]) != 1) {
      mexErrMsgTxt("Keywords must be strings.");
    }
    std::string keyword {mxArrayToString(prhs[i])};
    if (keyword == "Speed") {
      F = parseSpeedFunc(prhs[i + 1]);
    } else if (keyword == "Slowness") {
      F = parseSlownessFunc(prhs[i + 1]);
    } else if (keyword == "h") {
      h = parseH(prhs[i + 1]);
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
                      marcher_type & type, double & h, speed_func & F) {
  if (nrhs >= 2) {
    h = parseH(prhs[1]);
  }
  if (nrhs >= 3) {
    F = parseSpeedFunc(prhs[2]);
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
 *         (the uniform grid spacing), a function handle (the slowness
 *         function), a string (which method to use)
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[]) {
  /**
   * Ensure the correct number of outputs.
   */
  if (nlhs > 1) {
    mexErrMsgTxt("Too many outputs.");
  }

  /**
   * Setup arguments from passed parameters.
   */
  marcher_type type = marcher_type::basic;
  double h = 1;
  speed_func F = defaultSpeedFunc;
  if (nrhs >= 2 && mxIsClass(prhs[1], "char") && mxGetM(prhs[1]) == 1) {
    parseKeywordArguments(nlhs, plhs, nrhs, prhs, type, h, F);
  } else {
    parseRegularArguments(nlhs, plhs, nrhs, prhs, type, h, F);
  }

  /**
   * Setup input arguments on our end.
   */
  size_t M = mxGetM(prhs[0]);
  size_t N = mxGetN(prhs[0]);
  
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

  fmm_mex(out, in, M, N, h, F, type);

  if (F != defaultSpeedFunc) {
    cleanupUserFunc();
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
