# cython: embedsignature=True
# cython: language_level=3

import numpy as np

from enum import Enum

from libc.stdlib cimport malloc, free

class Neighborhood(Enum):
    OLIM4 = 0
    OLIM8 = 1
    OLIM6 = 2
    OLIM18 = 3
    OLIM26 = 4
    OLIM3D = 5

class Quadrature(Enum):
    MP0 = 0
    MP1 = 1
    RHR = 2

cdef extern from "olim_wrapper.h":
    enum neighborhood:
        OLIM4
        OLIM8
        OLIM6
        OLIM18
        OLIM26
        OLIM3D

    enum cost_func:
        MP0
        MP1
        RHR

    struct olim_wrapper_params:
        neighborhood nb
        cost_func F
        double h
        int *dims
        int ndims

    struct olim_wrapper:
        pass

    enum status:
        SUCCESS

    status olim_wrapper_init(olim_wrapper**, olim_wrapper_params*)
    status olim_wrapper_deinit(olim_wrapper**)
    status olim_wrapper_run(olim_wrapper*)
    status olim_wrapper_add_bd(olim_wrapper*, int*, double)
    status olim_wrapper_get_U_ptr(olim_wrapper*, double**)
    status olim_wrapper_get_s_ptr(olim_wrapper*, double**)
    status olim_wrapper_get_state_ptr(olim_wrapper*, char**)

cdef class Olim:
    cdef:
        olim_wrapper_params _p
        olim_wrapper* _w

    @property
    def dims(self):
        if self._p.ndims == 2:
            return self._p.dims[0], self._p.dims[1]
        elif self._p.ndims == 3:
            return self._p.dims[0], self._p.dims[1], self._p.dims[2]

    @property
    def U(self):
        cdef double * U
        cdef double[:, ::1] mv_2
        cdef double[:, :, ::1] mv_3
        olim_wrapper_get_U_ptr(self._w, &U)
        if self._p.ndims == 2:
            M, N = self.dims
            mv_2 = <double[:M, :N]> U
            return np.asarray(mv_2)
        elif self._p.ndims == 3:
            M, N, P = self.dims
            mv_3 = <double[:M, :N, :P]> U
            return np.asarray(mv_3)

    @U.setter
    def U(self, U_mv):
        if self._p.ndims == 2:
            self.set_U_2(U_mv)
        if self._p.ndims == 3:
            self.set_U_3(U_mv)

    cdef set_U_2(self, double[:, ::1] U_mv):
        cdef double * U = NULL
        olim_wrapper_get_U_ptr(self._w, &U)
        M, N = self.dims
        cdef double[:, ::1] mv = <double[:M, :N]> U
        mv[...] = U_mv

    cdef set_U_3(self, double[:, :, ::1] U_mv):
        cdef double * U = NULL
        olim_wrapper_get_U_ptr(self._w, &U)
        M, N, P = self.dims
        cdef double[:, :, ::1] mv = <double[:M, :N, :P]> U
        mv[...] = U_mv

    @property
    def s(self):
        cdef double * s
        cdef double[:, ::1] mv_2
        cdef double[:, :, ::1] mv_3
        olim_wrapper_get_s_ptr(self._w, &s)
        if self._p.ndims == 2:
            M, N = self.dims
            mv_2 = <double[:M, :N]> s
            return np.asarray(mv_2)
        elif self._p.ndims == 3:
            M, N, P = self.dims
            mv_3 = <double[:M, :N, :P]> s
            return np.asarray(mv_3)

    @s.setter
    def s(self, s_mv):
        if self._p.ndims == 2:
            self.set_s_2(s_mv)
        if self._p.ndims == 3:
            self.set_s_3(s_mv)

    cdef set_s_2(self, double[:, ::1] s_mv):
        cdef double * s = NULL
        olim_wrapper_get_s_ptr(self._w, &s)
        M, N = self.dims
        cdef double[:, ::1] mv = <double[:M, :N]> s
        mv[...] = s_mv

    cdef set_s_3(self, double[:, :, ::1] s_mv):
        cdef double * s = NULL
        olim_wrapper_get_s_ptr(self._w, &s)
        M, N, P = self.dims
        cdef double[:, :, ::1] mv = <double[:M, :N, :P]> s
        mv[...] = s_mv

    @property
    def state(self):
        cdef char * state
        cdef char[:, ::1] mv_2
        cdef char[:, :, ::1] mv_3
        olim_wrapper_get_state_ptr(self._w, &state)
        if self._p.ndims == 2:
            M, N = self.dims
            mv_2 = <char[:M, :N]> state
            return np.asarray(mv_2)
        elif self._p.ndims == 3:
            M, N, P = self.dims
            mv_3 = <char[:M, :N, :P]> state
            return np.asarray(mv_3)

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

    def __dealloc__(self):
        free(self._p.dims)
        err = olim_wrapper_deinit(&self._w)
        if err != SUCCESS:
            raise Exception('error!')

    def run(self):
        olim_wrapper_run(self._w)

    def add_bd(self, *args):
        if self._p.ndims == 2:
            if len(args) == 2:
                self.add_bd_2(list(args[:2]))
            elif len(args) == 3:
                self.add_bd_2(list(args[:2]), U=args[2])
            else:
                raise Exception('error!')
        elif self._p.ndims == 3:
            if len(args) == 3:
                self.add_bd_3(list(args[:3]))
            elif len(args) == 4:
                self.add_bd_3(list(args[:3]), U=args[3])
            else:
                raise Exception('error!')

    cdef add_bd_2(self, list inds, double U=0):
        cdef int inds_[2]
        inds_[0] = inds[0]
        inds_[1] = inds[1]
        olim_wrapper_add_bd(self._w, inds_, U)

    cdef add_bd_3(self, list inds, double U=0):
        cdef int inds_[3]
        inds_[0] = inds[0]
        inds_[1] = inds[1]
        inds_[2] = inds[2]
        olim_wrapper_add_bd(self._w, inds_, U)
