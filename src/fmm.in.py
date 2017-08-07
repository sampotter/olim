import numpy as _np

from cffi import FFI as _FFI

_ffi = _FFI()
_ffi.cdef("""
void fmm(double * out, bool * in, int M, int N, double h, double * S, 
         marcher_type type);
""")

_fmm = _ffi.dlopen("@LIBFMM_PATH@")

def _get_S_ptr(S, m, n):
    S_arr = np.ones((m, n))
    if S:
        for i in range(m):
            for j in range(n):
                S_arr[i, j] = S(i, j)
    return _ffi.from_buffer(S_arr)

def _get_marcher_type(method):
    return {
        'basic': 0, 'olim4_mp0': 1, 'olim4_rhr': 2, 'olim4_rhr_lut': 3,
        'olim8_mp0c': 4, 'olim8_mp1': 5, 'olim8_rhr': 6
    }[method]

def fmm(B, h=1, S=None, method='basic'):
    m, n = B.shape
    size = m*n
    in_ptr = _ffi.cast("bool *", _ffi.from_buffer(B))
    out_ptr = _ffi.new("double[%d]" % size)
    S_ptr = _get_S_ptr(S, m, n)
    marcher_type = _get_marcher_type(method)
    _fmm.fmm(out_ptr, in_ptr, m, n, h, marcher_type)
    return _np.frombuffer(_ffi.buffer(out_ptr, 8*size), dtype=_np.float64)
