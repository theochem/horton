cimport gbw

cdef extern from "horton/gbasis/cholesky.h":
    long cholesky(gbw.GB4IntegralWrapper* gbw4, double** uninit_result, double threshold)
