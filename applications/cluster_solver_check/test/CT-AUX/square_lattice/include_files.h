// File: include_files.h
//
// The value of QMC_INTEGRATOR_BIT is set by the cmake configuration.

// NFFT <--> DFT
const static bool DO_NFFT = true;

// random number generator
#include "dca/math_library/random_number_library/ranq2.hpp"
using random_number_generator = rng::ranq2;

// parallelization
#include "comp_library/parallelization_library/include_parallelization_library.h"

// include parameters
#include "phys_library/parameters/Parameters.h"

#include "phys_library/DCA+_data/DCA_data.h"

#include "phys_library/DCA+_loop/DCA+_loop_data.h"

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_cluster_solver.h"


#include "comp_library/IO_library/JSON/JSON_writer.h"