//-*-C++-*-


#include <cmath>
#include <cstdlib> 
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <complex>
#include <stdexcept>
#include <assert.h>
#include <limits>
#include <algorithm>
#include <typeinfo>
#include <pthread.h>
//#include <omp.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <queue>
#include <bitset>

using namespace std;

// typelist-operations

#include "type_list.h"   
#include "comp_library/type_list/type_list_definitions.h" 
using namespace TL;

// libraries

#include "fftw3.h"
#include "nfft3.h"

// enumerations

enum cluster_shape {BETT_CLUSTER, PARALLELEPIPED};
typedef cluster_shape cluster_shape_type;

enum mesh_shape {PARALLELLOGRAM, FIRST_BRILLOUIN_ZONE};
const static mesh_shape MESH_SHAPE = FIRST_BRILLOUIN_ZONE;//PARALLELLOGRAM; 

enum BRILLOUIN_ZONE {BRILLOUIN_ZONE_CUT_TEMPLATE,
		     FERMI_SURFACE_SQUARE_2D_LATTICE, SQUARE_2D_LATTICE,BODY_CENTERED_TETRAGONAL_A, BODY_CENTERED_TETRAGONAL_B, SIMPLE_TETRAGONAL, 
		     TRICLINIC, FACE_CENTERED_CUBIC, BODY_CENTERED_CUBIC, SIMPLE_CUBIC, HEXAGONAL, RHOMBOHEDRAL_A,
		     RHOMBOHEDRAL_B, SIMPLE_MONOCLINIC, ONE_FACE_CENTERED_MONOCLINIC_A, ONE_FACE_CENTERED_MONOCLINIC_B,
		     SIMPLE_ORTHOROMBIC, BASE_CENTERED_ORTHORHOMBIC, BODY_CENTERED_ORTHOROMBIC, ALL_FACE_CENTERED_ORTHORHOMBIC_A, ALL_FACE_CENTERED_ORTHORHOMBIC_B};

typedef BRILLOUIN_ZONE BRILLOUIN_ZONE_CUT_TYPE;


// provenance

//#include "provenance.h"

// include generic-methods
#include "include_generic_methods.h"

// blas/lapack
#include "include_blas_lapack_plans.h"
#include "include_linalg.h"

// various
#include "static_functions.h"
#include "include_symmetry_library.h"

// include function-library
#include "include_function_library.h"

// IO-library
#include "include_IO_operations.h"

// include plotting
#include "include_plotting.h"

// random number generator
#include "include_random_number_generator.h"

// parallelization
//#include "include_concurrency_serial.h"
//#include "include_concurrency_mpi.h"
//#include "include_concurrency_mpi_ft.h"
//#include "include_concurrency_ompi_ft.h"
#include "include_parallelization_library.h"

// profiling
#include "include_profiling.h"

#include "include_math_library.h"

// include domains
#include "include_DCA_domains.h"

// include VASP

#include "include_vasp.h"
