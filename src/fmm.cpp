#include "mex.h"

#include <iostream>
#include <string>
#include <unordered_map>

#include "fmm.h"

/**
 * These structs contain the parameters that need to be initialized in
 * order to instantiate one of the C++ fast marcher classes.
 */
struct parameters {
  marcher_type type;
  mxArray * S;
  double h, x0, y0;
  int M, N;
};
struct parameters_3d {
  parameters_3d(size_t const * dims): dims {dims} {}
  marcher_type type {BASIC};
  mxArray * S {nullptr};
  double h {1}, x0 {0}, y0 {0}, z0 {0};
  size_t const * dims;
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

static mxArray * getDefaultSMatrix3d(size_t const * dims) {
  mxArray * S = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  double * Spr = mxGetPr(S);
  int nelts = dims[0]*dims[1]*dims[2];
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
 * Fill the S matrix with the values of the user supplied speed
 * function.
 */
static void parseSpeedFunc3d(mxArray const * arg, parameters_3d & p) {
  if (!mxIsClass(arg, "function_handle")) {
    mexErrMsgTxt("Speed function argument must be a function.");
  }

  mxArray * X = mxCreateNumericArray(3, p.dims, mxDOUBLE_CLASS, mxREAL);
  mxArray * Y = mxCreateNumericArray(3, p.dims, mxDOUBLE_CLASS, mxREAL);
  mxArray * Z = mxCreateNumericArray(3, p.dims, mxDOUBLE_CLASS, mxREAL);

  double * Xpr = (double *) mxGetPr(X);
  double * Ypr = (double *) mxGetPr(Y);
  double * Zpr = (double *) mxGetPr(Z);

  int l = 0;
  for (int i = 0; i < p.dims[0]; ++i) {
    double y = p.h*i - p.y0;
    for (int j = 0; j < p.dims[1]; ++j) {
      double x = p.h*j - p.x0;
      for (int k = 0; k < p.dims[2]; ++k) {
        Xpr[l] = x;
        Ypr[l] = y;
        Zpr[l++] = p.h*k - p.z0;
      }
    }
  }

  mxArray * plhs[] = {nullptr};
  mxArray * prhs[] = {const_cast<mxArray *>(arg), X, Y, Z};
  mxArray * e = mexCallMATLABWithTrap(1, plhs, 4, prhs, "feval");
  if (e != nullptr) {
    std::cout << mxArrayToString(mxGetProperty(e, 0, "identifier")) << ": "
              << mxArrayToString(mxGetProperty(e, 0, "message")) << std::endl;
    std::exit(EXIT_FAILURE);
  }

  mxDestroyArray(X);
  mxDestroyArray(Y);
  mxDestroyArray(Z);
  mxDestroyArray(p.S); // free this before we lose our pointer to it

  p.S = plhs[0];
}

/**
 * Retrieve the marcher_type value from the MATLAB string passed as an
 * argument.
 */
static marcher_type parseMarcherType(mxArray const * arg) {
  auto const str = parseString(arg, "Method argument must be a string.");
  marcher_type type;
  if (str == "basic") type = BASIC;
  else if (str == "olim4_mp0") type = OLIM4_MP0;
  else if (str == "olim4_rhr") type = OLIM4_RHR;
  else if (str == "olim4_rhr_lut") type = OLIM4_RHR_LUT;
  else if (str == "olim8_mp0") type = OLIM8_MP0;
  else if (str == "olim8_mp1_bsearch") type = OLIM8_MP1_BSEARCH;
  else if (str == "olim8_mp1_gsl") type = OLIM8_MP1_GSL;
  else if (str == "olim8_mp1_secant") type = OLIM8_MP1_SECANT;
  else if (str == "olim8_rhr") type = OLIM8_RHR;
  else mexErrMsgTxt(("Invalid marcher type: " + str).c_str());
  return type;
}

/**
 * Iterate through passed keyword arguments and fill the parameter
 * struct.
 */
static void
parseKeywordArguments(mxArray * plhs[], int nrhs, mxArray const * prhs[],
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
parseRegularArguments(mxArray * plhs[], int nrhs, mxArray const * prhs[],
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

static void
parseAndRunMarcher(mxArray * plhs[], int nrhs, mxArray const * prhs[]) {
  /**
   * Get the height (M) and width (N) of the input boundary matrix.
   */
  int M = mxGetM(prhs[0]);
  int N = mxGetN(prhs[0]);

  /**
   * Set up arguments from passed parameters.
   */
  parameters p {BASIC, getDefaultSMatrix(M, N), 1, 0, 0, M, N}; // TODO: add ctor
  if (nrhs >= 2 && mxIsClass(prhs[1], "char") && mxGetM(prhs[1]) == 1) {
    parseKeywordArguments(plhs, nrhs, prhs, p);
  } else {
    parseRegularArguments(plhs, nrhs, prhs, p);
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

  fmm(out, in, M, N, p.h, mxGetPr(p.S), p.type);

  /**
   * Clean up.
   */
  free(in);
  mxDestroyArray(p.S);
}

/**
 * Fill 3D parameter struct from keyword arguments.
 */
static void
parseKeywordArguments3d(mxArray * plhs[], int nrhs, mxArray const * prhs[],
                        parameters_3d & p) {
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

  if (args.find("z0") != args.end())
    p.z0 = parseScalar(args["z0"], "Keyword argument 'z0' requires a scalar.");

  if (args.find("Speed") != args.end())
    parseSpeedFunc3d(args["Speed"], p);
}

/**
 * Iterate through fixed position arguments and fill the parameter
 * struct.
 */
static void
parseRegularArguments3d(mxArray * plhs[], int nrhs, mxArray const * prhs[],
                        parameters_3d & p) {
  if (nrhs >= 2) {
    p.h = parseScalar(
      prhs[1], "Second argument 'h' requires a scalar parameter.");
  }
  if (nrhs >= 3) {
    parseSpeedFunc3d(prhs[2], p);
  }
  if (nrhs >= 4) {
    p.type = parseMarcherType(prhs[3]);
  }
}

static void
parseAndRunMarcher3d(mxArray * plhs[], int nrhs, mxArray const * prhs[]) {
  /**
   * Get the dimensions of the boundary matrix.
   */
  size_t const * dims {mxGetDimensions(prhs[0])};
  size_t nelts {dims[0]*dims[1]*dims[2]};

  /**
   * Set up arguments from passed parameters.
   */
  parameters_3d p {dims};
  p.S = getDefaultSMatrix3d(dims);
  if (nrhs >= 2 && mxIsClass(prhs[1], "char") && mxGetM(prhs[1]) == 1) {
    parseKeywordArguments3d(plhs, nrhs, prhs, p);
  } else {
    parseRegularArguments3d(plhs, nrhs, prhs, p);
  }

  /**
   * Set up input arguments.
   */
  mxLogical * logicals = mxGetLogicals(prhs[0]);
  bool * in = (bool *) calloc(nelts, sizeof(bool));
  for (int i = 0; i < nelts; ++i) {
    in[i] = logicals[i];
  }

  /**
   * Initialize output and get pointer.
   */
  plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  double * out = static_cast<double *>(mxGetPr(plhs[0]));
 
  int dims_[] = {
    static_cast<int>(dims[0]),
    static_cast<int>(dims[1]),
    static_cast<int>(dims[2])
  };

  fmm3d(out, in, dims_, p.h, static_cast<double *>(mxGetPr(p.S)), p.type);

  /**
   * Clean things up.
   */
  free(in);
  mxDestroyArray(p.S);
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

  // TODO: do array error checking here (i.e. validate prhs[0])

  /**
   * Extract the number of dimensions of the MATLAB array to determine
   * whether to use a 2D or 3D marcher.
   */
  int ndims = static_cast<int>(mxGetNumberOfDimensions(prhs[0]));

  if (ndims == 2) {
    parseAndRunMarcher(plhs, nrhs, prhs);
  } else if (ndims == 3) {
    parseAndRunMarcher3d(plhs, nrhs, prhs);
  } else {
    mexErrMsgTxt("Boundary array must be either 2D or 3D.");
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
