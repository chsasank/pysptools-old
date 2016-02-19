#
#------------------------------------------------------------------------------
# Copyright (c) 2013-2014, Christian Therien
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#------------------------------------------------------------------------------
#
# nfindr.pyx - This file is part of the PySptools package.
#

"""
NFINDR function
"""


import numpy as np
cimport numpy as np
import scipy as sp

import math
import random
import pysptools.material_count.vd as vd
import eea

cimport cython

from libc.stdlib cimport malloc, free
from libc.float cimport FLT_MIN

ctypedef np.float64_t DTYPE_t
#ctypedef np.float32_t DTYPE_t

ctypedef DTYPE_t content
ctypedef content** mat

# LU decomposition implementation
# Pure C style
#
#

@cython.boundscheck(False)
@cython.wraparound(False)
cdef mat mat_new(int n) nogil:
    cdef int i
    cdef mat x = <mat>malloc(sizeof(content*) * n)

    for i in range(n):
      x[i] = <content*>malloc(n * sizeof(content));
    mat_zero(x, n);
    return x


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void mat_del(mat x, int n) nogil:
    for i in range(n):
        free(x[i])
    free(x)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void mat_zero(mat x, int n) nogil:
    cdef int i,j
    for i in range(n):
        for j in range(n):
            x[i][j] = 0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef mat mat_copy(void *s, int n) nogil:
    cdef int i,j
    cdef mat x = mat_new(n)
    cdef mat m = <mat>s
    for i in range(n):
        for j in range(n):
            x[i][j] = m[i][j]
    return x


@cython.boundscheck(False)
@cython.wraparound(False)
cdef mat_show(mat x, mid, n):
    print mid, '=',
    for i in range(n):
        if i == 0:
            print ' [ ',
        else:
            print "        ",
        for j in range(n):
            print x[i][j],
            if j < n - 1:
                print '  ',
            else:
                if i == n - 1:
                    print ' ]'
                else:
                    print


@cython.boundscheck(False)
@cython.wraparound(False)
cdef mat mat_mul(mat a, mat b, int n) nogil:
    cdef int i,j,k
    cdef mat c = mat_new(n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c[i][j] += a[i][k] * b[k][j]
    return c


@cython.boundscheck(False)
@cython.wraparound(False)
cdef content d_abs(content n) nogil:
    if n >= 0: return n
    return n * <content>-1.0


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void mat_pivot(mat a, mat p, int n) nogil:
    cdef int i,j,k,max_j
    cdef content tmp

    for i in range(n):
        for j in range(n):
            if i == j:
                p[i][j] = 1
            else:
                p[i][j] = 0
    for i in range(n):
        max_j = i
        for j in range(i,n):
            if (d_abs(a[j][i]) > d_abs(a[max_j][i])): max_j = j

        if (max_j != i):
            for k in range(n):
                tmp = p[i][k]
                p[i][k] = p[max_j][k]
                p[max_j][k] = tmp


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void mat_LU(mat A, mat L, mat U, mat P, int n) nogil:
    cdef int i,j,k
    cdef mat Aprime
    cdef content s

    mat_zero(L, n)
    mat_zero(U, n)
    mat_pivot(A, P, n)

    Aprime = mat_mul(P, A, n)

    for i in range(n): L[i][i] = 1

    for i in range(n):
        for j in range(n):
            if (j <= i):
                s = 0;
                for k in range(j): s+= L[j][k] * U[k][i]
                U[j][i] = Aprime[j][i] - s

            if (j >= i):
                s = 0
                for k in range(i) : s+= L[j][k] * U[k][i]
                if d_abs(U[i][i]) > FLT_MIN:
                    #L[j][i] = (Aprime[j][i] - s) / FLT_MIN
                    L[j][i] = (Aprime[j][i] - s) / U[i][i]

    mat_del(Aprime, n)

#
#
# end of LU decomposition


# C style OO interface to LU decomposition
#
#

cdef struct LU:
    int n
    mat A
    mat P
    mat L
    mat U


cdef LU * LU_new(int n) nogil:
    cdef LU *self = <LU*>malloc(sizeof(LU))
    self.P = mat_new(n)
    self.L = mat_new(n)
    self.U = mat_new(n)
    self.n = n

    return self


cdef void LU_decompose(LU *self, void *A) nogil:
    self.A = <mat>A
    mat_LU(self.A, self.L, self.U, self.P, self.n);


cdef mat LU_get_L(LU *self) nogil:
    return self.L


cdef mat LU_get_A(LU *self) nogil:
    return self.A


cdef mat LU_get_U(LU *self) nogil:
    return self.U


cdef void LU_show(LU *self, char x):
    if (x == 'A'):
        mat_show(self.A, 'A', self.n)
    if (x == 'P'):
        mat_show(self.P, 'P', self.n)
    if (x == 'L'):
        mat_show(self.L, 'L', self.n)
    if (x == 'U'):
        mat_show(self.U, 'U', self.n)


cdef void LU_free(LU *self) nogil:
    mat_del(self.P, self.n)
    mat_del(self.L, self.n)
    mat_del(self.U, self.n)
    free(self)

#
#
# end of C style OO interface


# Utility functions
#
#


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void mat_load(mat m, DTYPE_t [:,:] data, int n) nogil:
    for i in range(n):
        for j in range(n):
            m[i][j] = data[i,j]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void mat_load_column(Py_ssize_t rowstart, Py_ssize_t rowstop,
                          mat m, Py_ssize_t column_id, DTYPE_t [:] column) nogil:
    cdef Py_ssize_t j = 0
    for i in range(rowstart, rowstop):
        m[i][column_id] = column[j]
        j += 1

#
#
# end utility functions


# Partial determinant function
#
#

@cython.boundscheck(False)
@cython.wraparound(False)
cdef content abs_det(LU *lu, mat m, int n) nogil:
    cdef mat U
    cdef content det
    cdef int i

    LU_decompose(lu, m)
    U = LU_get_U(lu)
    det = 1
    for i in range(n):
        det = U[i][i] * det
    return d_abs(det)

#
#
# end determinant


@cython.boundscheck(False)
@cython.wraparound(False)
def NFINDR(data, p, transform=None, maxit=None, ATGP_init=False):
    print 'Use NFINDR.pyx'
    # data size
    nsamples, nvariables = data.shape

    if maxit == None:
        maxit = 3*p

    if transform == None:
        # transform as shape (N x p)
        transform = data
        transform = eea._PCA_transform(data, p-1)
    else:
        transform = transform
##    # data_pca as shape (N x p)
##    data_pca = eea._PCA_transform(data, p-1)
    transform = np.array(transform, copy=True, order='c', dtype=np.float)

    # Initialization
    # TestMatrix is a square matrix, the first row is set to 1
#    TestMatrix = np.zeros((p, p), dtype=np.float32)
    TestMatrix = np.zeros((p, p), dtype=np.float, order='c')
    TestMatrix[0,:] = 1
    IDX = None
    if ATGP_init == True:
        induced_em, idx = _ATGP(transform, p)
        IDX = np.array(idx, dtype=np.int64, order='c')
        for i in xrange(p):
            TestMatrix[1:p, i] = induced_em[i]
    else:
        IDX =  np.zeros((p), dtype=np.int64, order='c')
        for i in xrange(p):
            idx = int(math.floor(random.random()*nsamples))
            TestMatrix[1:p, i] = transform[idx]
            IDX[i] = idx

    IDX, it = _NFINDR(transform, TestMatrix, IDX, nsamples, p, maxit)

    E = np.zeros((len(IDX), nvariables), dtype=np.float, order='c')
    Et = np.zeros((len(IDX), p-1), dtype=np.float32)
    for j in xrange(len(IDX)):
        E[j] = data[IDX[j]]
        Et[j] = transform[IDX[j]]

    return E, Et, IDX, it


@cython.boundscheck(False)
@cython.wraparound(False)
def _NFINDR(np.ndarray[DTYPE_t, ndim=2, mode='c'] data_pca,
            np.ndarray[DTYPE_t, ndim=2, mode='c'] TestMatrix,
            np.ndarray[np.int64_t, ndim=1, mode='c'] IDX,
            Py_ssize_t nsamples,
            int p,
            int maxit):

    cdef int it
    cdef Py_ssize_t i, k
    cdef DTYPE_t actualVolume, volume, v1, v2
    cdef mat cTestMatrix

    cdef DTYPE_t [:,:] data_pca_view = data_pca
    cdef DTYPE_t [:,:] TestMatrix_view = TestMatrix
    cdef Py_ssize_t [:] IDX_view = IDX

    cdef LU* lu = LU_new(p)
    cTestMatrix = mat_new(p)
    mat_load(cTestMatrix, TestMatrix_view, p)

    v1 = -1
    v2 = 0
    it = 0
    actualVolume = 0
    while it < maxit and v2 > v1:
        for k in range(p):
            for i in range(nsamples):
                mat_load_column(1, p, cTestMatrix, k, data_pca_view[i])
                volume = abs_det(lu, cTestMatrix, p)
                if volume > actualVolume:
                    actualVolume = volume
                    IDX_view[k] = i
            mat_load_column(1, p, cTestMatrix, k, data_pca_view[IDX_view[k]])
        it += 1
        v1 = v2
        v2 = actualVolume

    LU_free(lu)
    mat_del(cTestMatrix, p)
    return IDX, it


@cython.boundscheck(False)
@cython.wraparound(False)
def _ATGP(np.ndarray[DTYPE_t, ndim=2, mode='c'] data, int p):
    cdef np.npy_intp nsamples = data.shape[0]
    cdef np.npy_intp nvariables = data.shape[1]
    cdef float max_energy, val
    cdef np.int64_t i
    cdef np.int64_t j
    cdef np.int64_t idx

    cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] E = np.zeros((p, nvariables), dtype=np.float)
    cdef np.ndarray[DTYPE_t, ndim=2, mode='c'] I = np.eye(nvariables, dtype=np.float)
    cdef np.ndarray[np.int64_t, ndim=1, mode='c'] IDX = np.zeros(p, dtype=np.int64)
    cdef np.ndarray[DTYPE_t, ndim=1, mode='c'] result = np.ndarray((nvariables,), dtype=np.float)

    # Algorithm initialization
    # the sample with max energy is selected as the initial endmember
    max_energy = -1
    idx = 0
    for j in range(nsamples):
        val = np.dot(data[j], data[j])
        if val > max_energy:
          max_energy = val
          idx = j

    # Initialization of the set of endmembers and the endmembers index vector
    E[0] = data[idx] # the first endmember selected
    IDX[0] = idx

    for i in range(p-1):
        UC = E[0:i+1]
        # Calculate the orthogonal projection with respect to the pixels at present chosen.
        # This part can be replaced with any other distance
        PU = I - np.dot(UC.T,np.dot(np.linalg.pinv(np.dot(UC,UC.T)),UC))
        max_energy = -1
        idx = 0
        # Calculate the most different pixel from the already selected ones according to
        # the orthogonal projection (or any other distance selected)
        for j in xrange(nsamples):
            r = data[j]
            np.dot(PU, r, out=result)
            val = np.dot(result.T, result)
            if val > max_energy:
                max_energy = val
                idx = j
    # The next chosen pixel is the most different from the already chosen ones
        E[i+1] = data[idx]
        # **fix2
        IDX[i+1] = idx
        #IDX[i] = idx

    return E, IDX

#-------------------------------------------------------------------------------
# Reference implementation that I use to compare to.
#
# The sFINDR function do an almost direct call to lapack.
# The sFINDR computing time is compared against the time used by the function
# NFINDR defined above in this file. The later is suppose to be faster. If not
# something is going wrong. By now, NFINDR is always faster than sFINDR when
# tested in the same context (same call to ATGP when it's time to take the initial
# working set).

@cython.boundscheck(False)
@cython.wraparound(False)
def sNFINDR(np.ndarray[np.float32_t, ndim=2] data_pca,
            np.ndarray[np.float32_t, ndim=2] TestMatrix,
            np.ndarray[np.int64_t, ndim=1] IDX,
            long int nsamples,
            int p,
            int maxit):

    cdef int k
    cdef Py_ssize_t i
    #cdef long int i
    cdef float actualVolume, volume, v2

    cdef Py_ssize_t it = 0
    cdef float v1 = -1.0

    actualVolume = 0
    v2 = actualVolume

    cdef float [:,:] data_pca_view = data_pca
    cdef float [:,:] TestMatrix_view = TestMatrix
    cdef long long [:] IDX_view = IDX

    while it < maxit and v2 > v1:
        for k in range(p):
            for i in range(nsamples):
                TestMatrix_view[1:p, k] = data_pca_view[i]
                volume = abs(sp.linalg._flinalg.sdet_r(TestMatrix)[0])
                if volume > actualVolume:
                    actualVolume = volume
                    IDX_view[k] = i
            TestMatrix_view[1:p, k] = data_pca_view[IDX_view[k]]
        it = it + 1
        v1 = v2
        v2 = actualVolume

    return IDX, it
