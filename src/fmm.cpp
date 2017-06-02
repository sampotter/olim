#include "mex.h"

#include <string>

#include "fmm_mex.hpp"

typedef double (*speed_func)(double, double);

double default_speed_func(double x, double y) {
	(void) x;
	(void) y;
	return 1;
}

static mxArray* user_speed_func_plhs[1] = {NULL};
static mxArray* user_speed_func_prhs[3] = {NULL, NULL, NULL};

double user_speed_func(double x, double y) {
	if (!mexCallMATLAB(1, user_speed_func_plhs, 3, user_speed_func_prhs,
					   "feval")) {
		// TODO: error handling
	}
	return mxGetScalar(user_speed_func_plhs[0]);
}

/**
 * nlhs -- number of outputs (expected: 1)
 * plhs -- outputs: a double matrix
 * nrhs -- number of inputs (expected: >=1)
 * prhs -- inputs: a logical matrix (the boundary points), a double
 *         (the uniform grid spacing), a function handle (the slowness
 *         function)
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, mxArray const * prhs[]) {
	/**
	 * Check to make sure the right number of arguments was passed.
	 */
	if (nlhs != 1) {
		mexErrMsgTxt("One output required.");
	}
	if (nrhs < 1) {
		mexErrMsgTxt("At least one input required.");
	}

	/**
	 * Make sure the second argument is a double scalar, if it was
	 * passed.
	 */
	if (nrhs >= 2 && (!mxIsClass(prhs[1], "double") || mxGetM(prhs[1]) > 1 ||
					  mxGetN(prhs[1]) > 1)) {
		mexErrMsgTxt("Second argument is not a double scalar.");
	}

	/**
	 * If passed, make sure the third argument is a function handle
	 * and do related MEX initialization for user_speed_func.
	 */
	if (nrhs >= 3) {
		if (!mxIsClass(prhs[2], "function_handle")) {
			mexErrMsgTxt("Third argument is not a function handle.");
		}

		user_speed_func_plhs[0] = mxCreateDoubleScalar(0);
		user_speed_func_prhs[0] = (mxArray *) prhs[2];
		user_speed_func_prhs[1] = mxCreateDoubleScalar(0);
		user_speed_func_prhs[2] = mxCreateDoubleScalar(0);
	}

	/**
	 * If it was passed, grab the string specifying the marcher to
	 * use.
	 */
	marcher_type type = marcher_type::basic;
	if (nrhs >= 4) {
		if (!mxIsClass(prhs[3], "char") || mxGetM(prhs[3]) != 1) {
			mexErrMsgTxt("Fourth argument must be a string.");
		}
		std::string str;
		str.reserve(mxGetN(prhs[3]));
		mxChar* chars = mxGetChars(prhs[3]);
		memcpy(&str[0], chars, mxGetN(prhs[3]));
		if (str == "basic") {
			type = marcher_type::basic;
		} else {
			mexErrMsgTxt("Invalid marcher type.");
		}
	}

	/**
	 * Setup input arguments on our end.
	 */
	size_t M = mxGetM(prhs[0]);
	size_t N = mxGetN(prhs[0]);
	mxLogical * logicals = mxGetLogicals(prhs[0]);
	bool * in = (bool *) calloc(M*N, sizeof(bool));
	for (int i = 0; i < M*N; ++i) in[i] = logicals[i];
	
	double h = nrhs >= 2 ? mxGetScalar(prhs[1]) : 1;
	speed_func F = nrhs >= 3 ? user_speed_func : default_speed_func;
	
	/**
	 * Initialize output and get pointer.
	 */
	plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
	double * out = mxGetPr(plhs[0]);

	fmm_mex(out, in, M, N, h, F, type);

	if (nrhs >= 3) {
		mxDestroyArray(user_speed_func_plhs[0]);
		mxDestroyArray(user_speed_func_prhs[1]);
		mxDestroyArray(user_speed_func_prhs[2]);
	}
}
