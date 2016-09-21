// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class manages the processor grouping for MPI.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_PROCESSOR_GROUPING_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_PROCESSOR_GROUPING_HPP

#include <cassert>
#include <mpi.h>

namespace dca {
namespace concurrency {
// dca::concurrency::

class MPIProcessorGrouping {
public:
  MPIProcessorGrouping() : id_(-1), nr_threads_(0) {}

  // We need a set-method since in ParallelizationMPI the constructor of this class is called before
  // MPI_Init.
  void set() {
    MPI_communication_ = MPI_COMM_WORLD;

    MPI_Comm_size(MPI_COMM_WORLD, &nr_threads_);
    MPI_Comm_rank(MPI_COMM_WORLD, &id_);
  }

  MPI_Comm get() const {
    return MPI_communication_;
  }
  int get_id() const {
    assert(id_ > -1);
    return id_;
  }
  int get_Nr_threads() const {
    assert(nr_threads_ > -1);
    return nr_threads_;
  }

  int first() const {
    return 0;
  }
  int last() const {
    return nr_threads_ - 1;
  }

private:
  MPI_Comm MPI_communication_;
  int id_;
  int nr_threads_;
};

}  // concurrency
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_PROCESSOR_GROUPING_HPP
