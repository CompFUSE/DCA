//-*-C++-*-

namespace COMP_LIB
{
  enum PARALLELIZATION_LIBRARY_NAMES {SERIAL_LIBRARY, POSIX_LIBRARY, OMP_LIBRARY, MPI_LIBRARY};
}

#include "type_map_interface.h"

#include "processor_grouping_interface.h"

#include "packing_interface.h"

#include "RNG_interface.h"

#include "collective_min_interface.h"
#include "collective_max_interface.h"
#include "collective_sum_interface.h"

#include "print_on_shell_interface.h"

#include "parallelization_template.h"

// MPI
#include <mpi.h>

#include "type_map_interface_mpi.h"

#include "processor_grouping_interface_mpi.h"

#include "packing_interface_mpi.h"
#include "collective_sum_interface_mpi.h"

#include "collective_min_interface_mpi.h"
#include "collective_max_interface_mpi.h"
#include "collective_sum_interface_mpi.h"

#include "parallelization_mpi.h"

// POSIX

#include <pthread.h>

#include "collective_sum_interface_pthreads.h"

#include "processor_grouping_interface_pthreads.h"

#include "parallelization_pthreads.h"

// thread-manager

#include "thread_manager_sum.h"
