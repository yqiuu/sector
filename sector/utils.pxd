from libc.stdlib cimport malloc, free
from libc.stdio cimport *
from libc.string cimport memcpy
from cextension cimport *


cdef int *init_1d_int(int[:] memview)
cdef float *init_1d_float(float[:] memview)
cdef double *init_1d_double(double[:] memview)
cdef llong_t *init_1d_llong(llong_t[:] memview)
