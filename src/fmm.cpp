#include "mex.h"

#include <iostream>
#include <string>
#include <unordered_map>

#include "fmm_mex.hpp"

/**
 * This struct contains the parameters that need to be initialized in
 * order to instantiate one of the C++ fast marcher classes.
 */
struct parameters {
  marcher_type type;
  mxArray * S;
  double h;
  double x0;
  double y0;
  int M;
  int N;
};

/**
 * The "S matrix" is a matrix of speed function values at each point
 * in the computational domain, *plus* the boundary values (which are
 * opaque to the user). This function fills this matrix with ones,
 * corresponding to the "default" speed function.
 */
static mxArray * getDefaultSMatrix(int M, int N) {
  mxArray * S = mxCreateDoubleMatrix(M + 2, N + 2, mxREAL);
  double * Spr = mxGetPr(S);
  int nelts = (M + 2)*(N + 2);
  for (int i = 0; i < nelts; ++i) {
    Spr[i] = 1;
  }
  return S;
}

/**
 * Retrieve a double scalar value from an mxArray, issue an error is
 * the mxArray isn't a double scalar.
 */
static double parseScalar(mxArray const * arg, const char * errMsgTxt) {
  if (!mxIsClass(arg, "double") || mxGetM(arg) > 1 || mxGetN(arg) > 1) {
    mexErrMsgTxt(errMsgTxt);
  }
  return mxGetScalar(arg);
}

/**
 * Retrieve a string from an mxArray, issue an error is the mxArray
 * isn't a MATLAB string.
 */
static std::string parseString(mxArray const * arg, const char * errMsgTxt) {
  if (!mxIsClass(arg, "char") || mxGetM(arg) != 1) {
    mexErrMsgTxt(errMsgTxt);
  }
  return mxArrayToString(arg);
}

/**
 * Fill the S matrix with the values of the user supplied speed
 * function.
 */
static void parseSpeedFunc(mxArray const * arg, parameters & p) {
  if (!mxIsClass(arg, "function_handle")) {
    mexErrMsgTxt("Speed function argument must be a function.");
  }

  mxArray * X = mxCreateDoubleMatrix(p.M, p.N, mxREAL);
  mxArray * Y = mxCreateDoubleMatrix(p.M, p.N, mxREAL);

  double * Xpr = mxGetPr(X), * Ypr = mxGetPr(Y);
  int k = 0;
  for (int i = 0; i < p.M; ++i) {
    double y = p.h*i - p.y0;
    for (int j = 0; j < p.N; ++j) {
      Xpr[k] = p.h*j - p.x0;
      Ypr[k++] = y;
    }
  }

  mxArray * plhs[] = {nullptr};
  mxArray * prhs[] = {const_cast<mxArray *>(arg), X, Y};
  mxArray * e = mexCallMATLABWithTrap(1, plhs, 3, prhs, "feval");
  if (e != nullptr) {
    std::cout << mxArrayToString(mxGetProperty(e, 0, "identifier")) << ": "
              << mxArrayToString(mxGetProperty(e, 0, "message")) << std::endl;
    std::exit(EXIT_FAILURE);
  }

  mxDestroyArray(X);
  mxDestroyArray(Y);
  mxDestroyArray(p.S); // free this before we lose our pointer to it

  p.S = plhs[0];
}

/**
 * Retrieve the marcher_type value from the MATLABstring passed as an
 * argument.
 */
static marcher_type parseMarcherType(mxArray const * arg) {
  auto const str = parseString(arg, "Method argument must be a string.");
  marcher_type type;
  if (str == "basic") type = marcher_type::basic;
  else if (str == "olim8_mp0") type = marcher_type::olim8_mp0;
  else if (str == "olim8_mp1") type = marcher_type::olim8_mp1;
  else if (str == "olim8_rhr") type = marcher_type::olim8_rhr;
  else mexErrMsgTxt(("Invalid marcher type: " + str).c_str());
  return type;
}

/**
 * Iterate through passed keyword arguments and fill the parameter
 * struct.
 */
static void
parseKeywordArguments(int nlhs, mxArray * plhs[],
                      int nrhs, mxArray const * prhs[],
                      parameters & p) {
  if (nrhs % 2 == 0) {
    mexErrMsgTxt("Keyword arguments must be passed in pairs.");
  }

  std::unordered_map<std::string, mxArray const *> args;
  for (int i = 1; i < nrhs; i += 2) {
    args[parseString(prhs[i], "Keywords must be strings.")] = prhs[i + 1];
  }

  if (args.find("Method") != args.end())
    p.type = parseMarcherType(args["Method"]);

  if (args.find("h") != args.end())
    p.h = parseScalar(args["h"], "Keyword argument 'h' requires a scalar.");

  if (args.find("x0") != args.end())
    p.x0 = parseScalar(args["x0"], "Keyword argument 'x0' requires a scalar.");

  if (args.find("y0") != args.end())
    p.y0 = parseScalar(args["y0"], "Keyword argument 'y0' requires a scalar.");

  if (args.find("Speed") != args.end())
    parseSpeedFunc(args["Speed"], p);
}

/**
 * Iterate through fixed position arguments and fill the parameter
 * struct.
 */
static void
parseRegularArguments(int nlhs, mxArray * plhs[],
                      int nrhs, mxArray const * prhs[],
                      parameters & p) {
  if (nrhs >= 2) {
    p.h = parseScalar(
      prhs[1], "Second argument 'h' requires a scalar parameter.");
  }
  if (nrhs >= 3) {
    parseSpeedFunc(prhs[2], p);
  }
  if (nrhs >= 4) {
    p.type = parseMarcherType(prhs[3]);
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
  parameters p {
    marcher_type::basic, getDefaultSMatrix(M, N), 1, 0, 0, M, N};
  if (nrhs >= 2 && mxIsClass(prhs[1], "char") && mxGetM(prhs[1]) == 1) {
    parseKeywordArguments(nlhs, plhs, nrhs, prhs, p);
  } else {
    parseRegularArguments(nlhs, plhs, nrhs, prhs, p);
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

  fmm_mex(out, in, M, N, p.h, mxGetPr(p.S), p.type);

  /**
   * Clean up.
   */
  free(in);
  mxDestroyArray(p.S);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
