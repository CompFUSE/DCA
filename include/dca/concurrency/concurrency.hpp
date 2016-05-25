// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file pulls in all concurrency header files.

#ifndef DCA_CONCURRENCY_CONCURRENCY_HPP
#define DCA_CONCURRENCY_CONCURRENCY_HPP

#include "parallelization_template.h"

// MPI
#ifdef MPI_SUPPORTED
#include "parallelization_mpi.h"
#endif

// POSIX threads
#include "parallelization_pthreads.h"

#include "thread_manager_sum.h"

#endif  // DCA_CONCURRENCY_CONCURRENCY_HPP
