// File: include_files.h

// Utilities library
#include <cstdlib>
#include <typeinfo>
#include <ctime>

// Numeric limits
#include <limits>

// Error handling
#include <stdexcept>
#include <cassert>

// Strings library
#include <cstring>

// Containers library
#include <vector>
#include <map>
#include <queue>

// Algorithms library
#include <algorithm>

// Numerics library
#include <cmath>
#include <complex>

// Input/output library
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>

// Parallelization libraries
#include <pthread.h>

// typelist-operations
#include "type_list.h"
using namespace TL;

// enumerations
#include "enumerations.hpp"

// BIT
const static bool QMC_INTEGRATOR_BIT = false;


// NFFT <--> DFT
const static bool DO_NFFT = true;

// blas/lapack
#include "comp_library/blas_lapack_plans/blas_lapack_plans.hpp"
#include "comp_library/linalg/linalg.hpp"

// various
#include "static_functions.h"
#include "include_symmetry_library.h"

// include generic-algorithms
#include "include_generic_methods.h"

// include function-library
#include "include_function_library.h"

// IO-library
#include "include_IO_operations.h"

// include plotting
#include "include_plotting.h"

// random number generator
#include "include_random_number_generator.h"

// include block-matrix operations

// parallelization
#include "include_parallelization_library.h"

// profiling
#include "include_profiling.h"

#include "include_math_library.h"

// include domains
#include "include_DCA_domains.h"

// type-dependent-conversions
#include "convert_DCA_types_to_index.h"

// include models
#include "include_Hamiltonians.h"

// include algorithms
#include "include_symmetries.h"
#include "compute_band_structure.h"

#include "include_DCA_steps.h"

// include parameters
#include "include_Parameters.h"

#include "include_DCA+_data.h"


#include "include_cluster_solver.h"

#include "include_DCA+_loop.h"

// analysis
#include "include_analysis.h"
