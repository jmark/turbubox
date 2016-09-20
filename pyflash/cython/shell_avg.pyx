cimport cython
cimport numpy as np
import numpy as np

from libc.math cimport sqrt

np.import_array()

@cython.boundscheck(False)
@cython.wraparound(False)
def shell_avg_3d(np.ndarray[double, ndim=3] X):
    cdef int Nx = X.shape[0]
    cdef int Ny = X.shape[1]
    cdef int Nz = X.shape[2]

    # cdef double[:] iws = np.zeros(int(np.sqrt(Nx*Nx/4+Ny*Ny/4+Nz*Nz)))
    # cdef double[:] nws = np.zeros(int(np.sqrt(Nx*Nx/4+Ny*Ny/4+Nz*Nz)))
    # cdef double[:] pws = np.zeros(int(np.sqrt(Nx*Nx/4+Ny*Ny/4+Nz*Nz)))

    cdef np.ndarray[double, ndim=1] iws = np.zeros(int(np.sqrt(Nx*Nx/4+Ny*Ny/4+Nz*Nz)))
    cdef np.ndarray[double, ndim=1] nws = np.zeros(int(np.sqrt(Nx*Nx/4+Ny*Ny/4+Nz*Nz)))
    cdef np.ndarray[double, ndim=1] pws = np.zeros(int(np.sqrt(Nx*Nx/4+Ny*Ny/4+Nz*Nz)))


    cdef int I,J,K,R
    cdef double r,n

    cdef int i,j,k

    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):

                I = i - Nx/2
                J = j - Ny/2
                K = k - Nz/2
                r = np.sqrt(I*I+J*J+K*K)
                R = int(r)

                iws[R] += r
                nws[R] += 1
                pws[R] += r*r * X[i,j,k]

    return np.asarray(iws),np.asarray(nws),np.asarray(pws)
