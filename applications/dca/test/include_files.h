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
// #include <omp.h>

// Deprecated libraries (<XXX.h> replaced by <cXXX>)
// #include <stdio.h>
// #include <math.h>
// #include <string.h>
// #include <stdlib.h>
// #include <assert.h>
// #include <time.h>

// Unused libaries
// #include <unistd.h>
// #include <dirent.h>

// #include <sys/resource.h>  // Placed into src/comp_library/profiler_library/events/time_file_name_changed.h
// #include <sys/time.h>      // Placed into performance_inspector_CPU.h and performance_inspector_GPU.h
// #include <bitset>          // Placed into advanced_fermionic_ed_type_definitions.h


// typelist-operations
#include <include/dca/util/type_list.hpp>
#include <include/dca/util/type_utils.hpp>
using namespace dca::util;

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
#include "math_library/static_functions.h"
#include "phys_library/domains/cluster/symmetries/include_symmetry_library.h"
//#include "statistical_methods.h"

// include function-library
#include "comp_library/function_library/include_function_library.h"

// IO-library
#include "comp_library/IO_library/include_IO_operations.h"

// include plotting
#include "comp_library/function_plotting/include_plotting.h"

// random number generator
#include "math_library/random_number_library/include_random_number_generator.h"

// include block-matrix operations
//#include "include_blocked_blas_calls.h"

// parallelization
#include "comp_library/parallelization_library/include_parallelization_library.h"

// profiling
#include "comp_library/profiler_library/include_profiling.h"

#include "math_library/math_library.hpp"

// include domains
#include "phys_library/domains/include_DCA_domains.h"
// #include "numerical_error_domain.h"
// #include "DCA_iteration_domain.h"
// #include "Feynman_expansion_order_domain.h"
// #include "electron_band_domain.h"
// #include "electron_spin_domain.h"
// #include "HS_spin_domain.h"
// #include "HS_field_sign_domain.h"
// #include "particle_number_domain.h"
// #include "brillouin_zone_cut_domain.h"
// #include "Brillouin_zone_cut.h"
// #include "legendre_domain.h"
// #include "point_group_operation_dmn.h"
// #include "time_domain.h"
// #include "time_domain_left_oriented.h"
// #include "frequency_domain.h"
// #include "frequency_domain_real_axis.h"
// #include "frequency_domain_imag_axis.h"
// #include "frequency_domain_compact.h"

// type-dependent-conversions
#include "phys_library/domains/convert_DCA_types_to_index.h"

// include models
#include "phys_library/parameters/models/include_Hamiltonians.h"
// #include "dft_model.h"
// #include "Koshevnikov_model.h"

// include algorithms
#include "phys_library/DCA+_step/symmetrization/include_symmetries.h"
// #include "include_tetrahedron_mesh.h"
// #include "include_symmetries.h"
#include "phys_library/DCA+_algorithms/compute_band_structure/compute_band_structure.h"
// #include "adjust_chemical_potential.h"

#include "phys_library/DCA+_step/include_DCA_steps.h"

// include parameters
#include "phys_library/parameters/include_Parameters.h"

#include "phys_library/DCA+_data/DCA_data.h"

// #include "include_DCA_steps.h"

#include "phys_library/DCA+_step/cluster_solver/include_cluster_solver.h"

#include "phys_library/DCA+_loop/DCA+_loop.h"

// include CPE
//#include "include_CPE.h"

// include ED
//#include "include_ed.h"

// include SE
//#include "include_series_expansion.h"

// analysis
#include "phys_library/DCA+_analysis/include_analysis.h"
