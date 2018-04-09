// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class manages the processor grouping for MPI.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_PROCESSOR_GROUPING_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_PROCESSOR_GROUPING_HPP

#include <cassert>
#include <vector>

#include <mpi.h>


namespace dca {
namespace parallel {
// dca::parallel::

class MPIProcessorGrouping {
public:
  MPIProcessorGrouping();

  ~MPIProcessorGrouping();


  inline MPI_Comm get() const {
    return MPI_communication_;
  }
  inline int get_id() const {
    assert(id_ > -1);
    return id_;
  }
  inline int get_world_id() const {
    assert(world_id_ > -1);
    return world_id_;
  }
  inline int get_Nr_threads() const {
    assert(nr_threads_ > -1);
    return nr_threads_;
  }

  inline int first() const {
    return 0;
  }
  inline int last() const {
    return nr_threads_ - 1;
  }
  inline bool isValid() const {
    return id_ >= 0;
  }

private:

  bool testValidity() const;

  void printRemovedProcesses(const std::vector<int>& valid_ids) const;

private:
  MPI_Group MPI_group_ = nullptr;
  MPI_Comm MPI_communication_ = nullptr;
  int world_id_ = -1;
  int id_ = -1;
  int nr_threads_ = -1;
  int world_size_ = -1;
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_PROCESSOR_GROUPING_HPP
