//-*-C++-*-                                                                                                                                                                                                                                                                         

#ifndef QMC_MATRIX_SCALARTYPE_GPU_H
#define QMC_MATRIX_SCALARTYPE_GPU_H

namespace LIN_ALG {

    template<>
    struct MATRIX_SCALARTYPE<float, GPU>
    {
	typedef float new_scalartype;
    };
    
    template<>
    struct MATRIX_SCALARTYPE<double, GPU>
    {
	typedef double new_scalartype;
    };
    
    template<>
    struct MATRIX_SCALARTYPE<std::complex<float>, GPU>
    {
	typedef cuComplex new_scalartype;
    };
    
    template<>
    struct MATRIX_SCALARTYPE<std::complex<double>, GPU>
    {
	typedef cuDoubleComplex new_scalartype;
    };
}

#endif
