cimport gbasis
cimport ints

cdef extern from "horton/gbasis/gbw.h":
    cdef cppclass GB4IntegralWrapper:
        GB4IntegralWrapper(gbasis.GOBasis* gobasis, ints.GB4Integral* gb4int)
        long get_nbasis()
        void select_2index(long index0, long index2,
                            long* pbegin0, long* pend0,
                            long* pbegin2, long* pend2)
        void compute()
        void compute_diagonal(double* diagonal)
        double* get_2index_slice(long index0, long index2)
