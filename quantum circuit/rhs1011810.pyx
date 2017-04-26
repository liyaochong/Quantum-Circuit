# This file is generated automatically by QuTiP.
# (C) 2011 and later, P. D. Nation & J. R. Johansson

import numpy as np
cimport numpy as np
cimport cython
from qutip.cy.spmatfuncs cimport spmvpy
from qutip.cy.interpolate cimport interp, zinterp
from qutip.cy.math cimport erf
cdef double pi = 3.14159265358979323

include '/home/chen/anaconda3/lib/python3.6/site-packages/qutip/cy/complex_math.pxi'

ctypedef np.complex128_t CTYPE_t
ctypedef np.float64_t DTYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
def cy_td_ode_rhs(
        double t,
        np.ndarray[CTYPE_t, ndim=1] vec,
        np.ndarray[CTYPE_t, ndim=1] data0,np.ndarray[int, ndim=1] idx0,np.ndarray[int, ndim=1] ptr0,
        np.ndarray[CTYPE_t, ndim=1] data1,np.ndarray[int, ndim=1] idx1,np.ndarray[int, ndim=1] ptr1,
        np.ndarray[CTYPE_t, ndim=1] data2,np.ndarray[int, ndim=1] idx2,np.ndarray[int, ndim=1] ptr2,
        np.ndarray[CTYPE_t, ndim=1] data3,np.ndarray[int, ndim=1] idx3,np.ndarray[int, ndim=1] ptr3,
        np.ndarray[CTYPE_t, ndim=1] data4,np.ndarray[int, ndim=1] idx4,np.ndarray[int, ndim=1] ptr4,
        float f0,
        float Omega0,
        int width0,
        float f1,
        float Omega1,
        int width1):
    
    cdef Py_ssize_t row
    cdef int num_rows = len(vec)
    cdef np.ndarray[CTYPE_t, ndim=1] out = np.zeros((num_rows),dtype=np.complex)
     
    spmvpy(data0, idx0, ptr0, vec, 1.0, out)
    spmvpy(data1, idx1, ptr1, vec, (0)*(0<t<=80), out)
    spmvpy(data2, idx2, ptr2, vec, (0)*(0<t<=80), out)
    spmvpy(data3, idx3, ptr3, vec, (-Omega0/2*(np.exp(-(t-20-0)**2/2.0/width0**2)*np.cos(t*f0+np.pi/2)-(t-20-0)/2.0/width0**2/-1.57079632679*np.exp(-(t-20-0)**2/2.0/width0**2)*np.cos(t*f0+np.pi)))*(0<t<=40)(Omega0*(np.exp(-(t-20-0-40)**2/2.0/width0**2)*np.cos(t*f0)-(t-20-0-40)/2.0/width0**2/-1.57079632679*np.exp(-(t-20-0-40)**2/2.0/width0**2)*np.cos(t*f0+np.pi/2.0)))*(40<t<=80), out)
    spmvpy(data4, idx4, ptr4, vec, (-Omega1/2*(np.exp(-(t-20-0)**2/2.0/width1**2)*np.cos(t*f1+np.pi/2)-(t-20-0)/2.0/width1**2/-1.57079632679*np.exp(-(t-20-0)**2/2.0/width1**2)*np.cos(t*f1+np.pi)))*(0<t<=40)(Omega1*(np.exp(-(t-20-0-40)**2/2.0/width1**2)*np.cos(t*f1)-(t-20-0-40)/2.0/width1**2/-1.57079632679*np.exp(-(t-20-0-40)**2/2.0/width1**2)*np.cos(t*f1+np.pi/2.0)))*(40<t<=80), out)
    return out
