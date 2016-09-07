// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This header file defines the supported underlying parallelization libraries.
// TODO: - Remove Pthreads and OpenMP?
//       - Make conform with coding style.

#ifndef DCA_CONCURRENCY_CONCURRENCY_TYPES_HPP
#define DCA_CONCURRENCY_CONCURRENCY_TYPES_HPP

namespace dca {
namespace concurrency {
// dca::concurrency::

enum PARALLELIZATION_LIBRARY_NAMES { SERIAL_LIBRARY, POSIX_LIBRARY, OMP_LIBRARY, MPI_LIBRARY };

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_CONCURRENCY_TYPES_HPP
