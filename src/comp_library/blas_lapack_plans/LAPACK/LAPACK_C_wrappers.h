
//-*-C++-*-
// ****************************************************************************
// * C++ wrapper for LAPACK                                                   *
// *                                                                          *
// * Thomas Schulthess, ORNL, October 1999                                    *
// * Peter Staar, ETHZ, December 2011                                         *
// ****************************************************************************

#ifndef LAPACK_C_WRAPPERS
#define LAPACK_C_WRAPPERS

/** \file LAPACK
 *  \author Thomas C. Schulthess, Michael S. Summers, Peter Staar
 */

#include <complex>

namespace LAPACK {

// ============================================================================


extern "C" float slamch_(char* JOBZ);
extern "C" double dlamch_(char* JOBZ);
// solve
extern "C" void sgesv_(int*, int*, float*, int*, int*, float*, int*, int*);
extern "C" void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
extern "C" void cgesv_(int*, int*, std::complex<float>*, int*, int*, std::complex<float>*, int*,
                       int*);
extern "C" void zgesv_(int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*,
                       int*);

} /* namespace LAPACK */

#endif  // PSIMAG_LAPACK_H
