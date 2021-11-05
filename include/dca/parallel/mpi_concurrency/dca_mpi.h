#include <mpi.h>

#ifdef MPICH_NUMVERSION
// Maybe this is fixed in mpich 4.0? Not checked yet.
#if MPICH_NUMVERSION < 40000000

/* Fix broken MPI-3 C++ types in mpich */
// Note: pre-processors only compare numbers, so the comparison
// has to be done as C/C++ code.
#undef MPI_CXX_BOOL
#define MPI_CXX_BOOL                MPI_C_BOOL

#undef MPI_CXX_FLOAT_COMPLEX
#define MPI_CXX_FLOAT_COMPLEX       MPI_C_FLOAT_COMPLEX

#undef MPI_CXX_DOUBLE_COMPLEX
#define MPI_CXX_DOUBLE_COMPLEX      MPI_C_DOUBLE_COMPLEX

#undef MPI_CXX_LONG_DOUBLE_COMPLEX
#define MPI_CXX_LONG_DOUBLE_COMPLEX MPI_C_LONG_DOUBLE_COMPLEX

#else
#warning "MPICH may have broken MPI_CXX_FLOAT_COMPLEX definitions."
#endif // MPICH_NUMVERSION
#endif // MPICH_NUMVERSION
