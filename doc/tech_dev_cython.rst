..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2016 The HORTON Development Team
    :
    : This file is part of HORTON.
    :
    : HORTON is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : HORTON is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --

Accelerating HORTON code with Cython
####################################

HORTON was designed to prioritize ease of programming over performance. This is a reasonable decision in light of the fact that the vast majority of time in a quantum chemistry code is only spent in a relatively small section of code. In HORTON, we rewrite these critical portion s of code in C++ and link them into Python using the Cython framework. We identify these critical portions using profiling tools.


Before you begin
================

There are several downsides to accelerating code with Cython. Please make sure they are acceptable to you before starting to code.

- Developer tools will break. PDB cannot read C++ debug symbols. cProfile will break. IDEs will not syntax highlight correctly. Valgrind will report false positive memory leaks.
- Cython is still a fairly new project. The API may not be stable. Don't be surprised if your code breaks after a few versions.
- A steep learning curve. Unless you are already familiar with C/C++ and Python profiling tools, you may not obtain the speed up you expected. Memory allocation and management of arrays is particularly tricky.


Basic example
=============

Background
----------

We will take a simplified example from the slicing code of the matrix class. The original Cholesky decomposed 2-electron integrals had an operation to slice along 3-indices which was consuming a significant portion of the time in the code. This was implemented using the ``numpy.einsum`` method.

The original code is here

.. code-block:: python

    def get_slice_slow(self):
        return numpy.einsum("ack, bck-> abc", B, B_prime)

This code takes a slice of B where the indices ``c`` are kept the same and then contracts across the last index.

.. math::
    A_{abc} = \sum_k B_{ack} B'_{bck}

A quick check using the python cProfile module (``python -m cProfile -o slice.pstats get_slice_slow.py; gprof2dot -f pstats slice.pstats | dot -Tpng -o slice.png``) showed ``get_slice_slow`` method required almost 40% of the total code runtime. Since this operation was simple to implement in C++, it was a good candidate for Cythonizing.

The C++ code to implement the same operation is below:

.. code-block:: c++

    //get_slice_abcc.cpp
    void get_slice_abcc(double* B, double* B_prime, double* out, long nbasis, long nvec){
        new long k; //example
        for (k=0; k<nvec; k++){
            for (long a=0; a<nbasis; a++){
                for (long b=0; b<nbasis; b++){
                    for (long c=0; c<nbasis; c++){
                        out[a*nbasis*nbasis + b*nbasis + c] += inp[k*nbasis*nbasis + a*nbasis + c] * inp2[k*nbasis*nbasis + b*nbasis + c];
                    }
                }
            }
        }
        delete k; //example
    }

and the header is below:

.. code-block:: c++

    //get_slice_abcc.h
    void get_slice_abcc(double* inp, double* inp2, double* out, long nbasis, long nvec);

This code now needs to be interfaced with Python using Cython.


Cythonizing your code
---------------------

First create a Cython `.pxd` header. This file provides information for Cython to link your compiled code to the cython file later.

.. code-block:: python

    #get_slice_abcc.pxd
    cdef extern from "get_slice_abcc.h":
        void get_slice_abcc(double* B, double* B_prime, double* out, long nbasis, long nvec)

You'll notice here that the Cython header is remarkably similar to the C++ header. There are a few keywords introduced here, the most significant being ``cdef``. The `.pxd` files are python syntax with a few other keywords and syntax for pointers. See the Cython documentation on C++ below for more details on how to Cythonize classes and more.

The `.pyx` file is where brunt of the work by Cython is done. It is also python syntax with a few extra keywords.

.. code-block:: python

    #cext.pyx
    cimport get_slice_abcc
    def get_slice_fast(np.ndarray[double, ndim=3] B not None,
                        np.ndarray[double, ndim=3] B_prime not None,
                        np.ndarray[double, ndim=3] out not None,
                       long nbasis, long nvec):

    assert B.flags['C_CONTIGUOUS']
    assert B.shape[0] == nvec
    assert B.shape[1] == B.shape[2] == nbasis
    #etc...

    get_slice_abcc.get_slice_abcc(&B[0, 0, 0], &B_prime[0, 0, 0], &out[0, 0, 0], nbasis, nvec)
    return out

There are several things to note here:

- The arguments are statically typed.
- The Numpy arrays have their datatypes declared as well as the number of dimensions
- It is good practice to have safety checks because the code in `.pyx` files will *not* give clean stack traces.
- Python and Numpy use ``long`` datatypes by default.
- You can pass the address of the first element of a Numpy array to a function expecting ``double*`` as long as it is contiguous.

There are several other nuances not illustrated in this example, but they are well covered in the Cython documentation below. Users should be particularly cognizant of whether variables are Python-style (dynamic typed) or C-style (static typed). In our example above, everything is static typed as the method declaration declares everything.


Additional notes
================

The above example leaves all memory management to the Python interpreter. This is not always possible, especially when implementing iterative algorithms in C/C++ code. There is no issue when memory is allocated and deallocated dynamically in the C++ code as in the example above. However, if memory must be allocated by C++ and freed by Python, it can be much more complicated. The reverse case, memory allocated by Python and freed by C++, should be much more rare and won't be covered here.

The most common form of memory allocated in C++ and passed back to Python for management is likely Numpy arrays. We will show a code snippet for managing this.

.. code-block:: python

    cdef double* data = NULL
    cdef np.npy_intp dims[3]

    nvec = calculate_cholesky(&data)
    dims[0] = <np.npy_intp> nvec
    dims[1] = <np.npy_intp> nbasis
    dims[2] = <np.npy_intp> nbasis

    result = numpy.PyArray_SimpleNewFromData(3, dims, np.NPY_DOUBLE, data)

The method ``PyArray_SimpleNewFromData`` creates a new Numpy array from memory which has already been allocated. The numpy data types must be specified, as well as the dimensionality. Data is simply a 1D ``double*`` array of size nvec * nbasis * nbasis.


Further reading
===============

http://docs.cython.org/

http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html
