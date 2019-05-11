# distutils: language=c++
# cython: embedsignature=True
# cython: language_level=3

import numpy as np

from enum import Enum

from libc.stdlib cimport malloc, free
from libcpp cimport bool

class Neighborhood(Enum):
    OLIM4 = 0
    OLIM8 = 1
    OLIM6 = 2
    OLIM18 = 3
    OLIM26 = 4
    OLIM3D = 5
    FMM2 = 6
    FMM3 = 7

class Quadrature(Enum):
    MP0 = 0
    MP1 = 1
    RHR = 2

class State(Enum):
    VALID = 0
    TRIAL = 1
    FAR = 2
    BOUNDARY = 3
    FREE = 4

cdef extern from "olim_wrapper.h":
    enum neighborhood:
        OLIM4
        OLIM8
        OLIM6
        OLIM18
        OLIM26
        OLIM3D
        FMM2
        FMM3

    enum cost_func:
        MP0
        MP1
        RHR

    enum status:
        SUCCESS

    struct fac_src_wrapper_params:
        int ndims
        double *coords
        double s

    struct fac_src_wrapper:
        pass

    status fac_src_wrapper_init(fac_src_wrapper**, fac_src_wrapper_params *);
    status fac_src_wrapper_deinit(fac_src_wrapper**);

    struct olim_wrapper_params:
        neighborhood nb
        cost_func F
        double h
        int *dims
        int ndims

    struct olim_wrapper:
        pass

    status olim_wrapper_init(olim_wrapper**, olim_wrapper_params*)
    status olim_wrapper_deinit(olim_wrapper**)
    status olim_wrapper_solve(olim_wrapper*)
    status olim_wrapper_step(olim_wrapper*, int*)
    status olim_wrapper_peek(olim_wrapper*, double*, int*, bool*)
    status olim_wrapper_adjust(olim_wrapper*, int*, double)
    status olim_wrapper_add_src(olim_wrapper*, int*, double)
    status olim_wrapper_add_bd(olim_wrapper*, int*)
    status olim_wrapper_add_free(olim_wrapper*, int*)
    status olim_wrapper_set_fac_src(olim_wrapper*, int*, void*)
    status olim_wrapper_get_U_ptr(olim_wrapper*, double**)
    status olim_wrapper_get_s_ptr(olim_wrapper*, double**)
    status olim_wrapper_get_state_ptr(olim_wrapper*, char**)

cdef class FacSrc:
    cdef:
        fac_src_wrapper_params _p
        fac_src_wrapper* _w

    def __cinit__(self, coords, s):
        self._p.ndims = len(coords)
        self._p.coords = <double*>malloc(self._p.ndims*sizeof(double))
        self._p.s = s
        err = fac_src_wrapper_init(&self._w, &self._p)
        if err != SUCCESS:
            raise Exception('error!')

    def __dealloc__(self):
        free(self._p.coords)
        fac_src_wrapper_deinit(&self._w)

cdef class Olim:
    cdef:
        olim_wrapper_params _p
        olim_wrapper* _w

    @property
    def quad(self):
        return Quadrature(self._p.F)

    @property
    def nb(self):
        return Neighborhood(self._p.nb)

    @property
    def h(self):
        return self._p.h

    @property
    def ndims(self):
        return self._p.ndims

    @property
    def shape(self):
        if self.ndims == 2:
            return self._p.dims[0], self._p.dims[1]
        elif self.ndims == 3:
            return self._p.dims[0], self._p.dims[1], self._p.dims[2]

    @property
    def size(self):
        return np.array(self.shape).prod()

    @property
    def U(self):
        cdef double* U
        olim_wrapper_get_U_ptr(self._w, &U)
        dims_ = np.array(self.shape) + 2*np.ones(self.ndims, dtype=np.int)
        size_ = dims_.prod()
        arr = np.asarray(<double[:size_]> U).reshape(dims_)
        if self.ndims == 2:
            return arr[1:-1, 1:-1]
        elif self.ndims == 3:
            return arr[1:-1, 1:-1, 1:-1]

    @U.setter
    def U(self, U_mv):
        if self._p.ndims == 2:
            self.set_U_2(U_mv)
        if self._p.ndims == 3:
            self.set_U_3(U_mv)

    cdef set_U_2(self, double[:, ::1] U_mv):
        cdef double * U = NULL
        olim_wrapper_get_U_ptr(self._w, &U)
        M, N = self.shape
        cdef double[:, ::1] mv = <double[:(M + 2), :(N + 2)]> U
        mv[1:-1, 1:-1] = U_mv

    cdef set_U_3(self, double[:, :, ::1] U_mv):
        cdef double * U = NULL
        olim_wrapper_get_U_ptr(self._w, &U)
        M, N, P = self.shape
        cdef double[:, :, ::1] mv = <double[:(M + 2), :(N + 2), :(P + 2)]> U
        mv[1:-1, 1:-1, 1:-1] = U_mv

    @property
    def s(self):
        cdef double* s
        olim_wrapper_get_s_ptr(self._w, &s)
        dims_ = np.array(self.shape) + 2*np.ones(self.ndims, dtype=np.int)
        size_ = dims_.prod()
        arr = np.asarray(<double[:size_]> s).reshape(dims_)
        if self.ndims == 2:
            return arr[1:-1, 1:-1]
        elif self.ndims == 3:
            return arr[1:-1, 1:-1, 1:-1]

    @s.setter
    def s(self, s_mv):
        if self._p.ndims == 2:
            self.set_s_2(s_mv)
        if self._p.ndims == 3:
            self.set_s_3(s_mv)

    cdef set_s_2(self, double[:, ::1] s_mv):
        cdef double * s = NULL
        olim_wrapper_get_s_ptr(self._w, &s)
        M, N = self.shape
        cdef double[:, ::1] mv = <double[:(M + 2), :(N + 2)]> s
        mv[1:-1, 1:-1] = s_mv

    cdef set_s_3(self, double[:, :, ::1] s_mv):
        cdef double * s = NULL
        olim_wrapper_get_s_ptr(self._w, &s)
        M, N, P = self.shape
        cdef double[:, :, ::1] mv = <double[:(M + 2), :(N + 2), :(P + 2)]> s
        mv[1:-1, 1:-1, 1:-1] = s_mv

    @property
    def state(self):
        cdef char * state
        olim_wrapper_get_state_ptr(self._w, &state)
        dims_ = np.array(self.shape) + 2*np.ones(self.ndims, dtype=np.int)
        size_ = dims_.prod()
        arr = np.asarray(<char[:size_]> state).reshape(dims_)
        if self.ndims == 2:
            return arr[1:-1, 1:-1]
        elif self.ndims == 3:
            return arr[1:-1, 1:-1, 1:-1]

    def __cinit__(self, nb, quad, s, double h):
        self._p.nb = nb.value
        self._p.F = quad.value
        self._p.h = h
        self._p.dims = <int *> malloc(s.ndim*sizeof(int))
        for i in range(s.ndim):
            self._p.dims[i] = s.shape[i]
        self._p.ndims = s.ndim

        err = olim_wrapper_init(&self._w, &self._p)
        if err != SUCCESS:
            raise Exception('error!')

        self.s[...] = s

# TODO: this causes some problems when an instance of Olim is managed
# as instance member of another class...? Not 100% sure...

#    def __dealloc__(self):
#        free(self._p.dims)
#        err = olim_wrapper_deinit(&self._w)
#        if err != SUCCESS:
#            raise Exception('error!')

    def solve(self):
        olim_wrapper_solve(self._w)

    def step(self):
        cdef int lin
        olim_wrapper_step(self._w, &lin)
        if lin == -1:
            return None
        else:
            return np.unravel_index(lin, self.shape)

    def min(self):
        cdef double value
        cdef bool empty
        cdef int lin
        olim_wrapper_peek(self._w, &value, &lin, &empty)
        if not empty:
            return value, np.unravel_index(lin, self.shape)

    # TODO: this can be simplified
    cdef get_inds_mv(self, inds):
        cdef int[::1] mv = np.empty((self._p.ndims,), dtype=np.intc)
        cdef int i
        for i, ind in enumerate(inds):
            mv[i] = ind
        return mv

    cpdef adjust(self, inds, U):
        cdef int[::1] mv = self.get_inds_mv(inds)
        olim_wrapper_adjust(self._w, &mv[0], U)

    cpdef add_src(self, inds, U=0):
        cdef int[::1] mv = self.get_inds_mv(inds)
        olim_wrapper_add_src(self._w, &mv[0], U)

    cpdef add_bd(self, inds):
        cdef int[::1] mv = self.get_inds_mv(inds)
        olim_wrapper_add_bd(self._w, &mv[0])

    cpdef add_free(self, inds):
        cdef int[::1] mv = self.get_inds_mv(inds)
        olim_wrapper_add_free(self._w, &mv[0])

    cpdef set_fac_src(self, inds, FacSrc fs):
        cdef int[::1] mv = self.get_inds_mv(inds)
        olim_wrapper_set_fac_src(self._w, &mv[0], fs._w)
