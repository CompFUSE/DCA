// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

  MPIProcessorGrouping(bool (*test)());

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
  inline int get_size() const {
    assert(size_ > -1);
    return size_;
  }

  inline int first() const {
    return 0;
  }
  inline int last() const {
    return size_ - 1;
  }
  inline bool isValid() const {
    return id_ >= 0;
  }

private:

  static bool defaultTest();

  void printRemovedProcesses(const std::vector<int>& valid_ids) const;

private:
  MPI_Group MPI_group_ = 0;
  MPI_Comm MPI_communication_ = 0;
  int world_id_ = -1;
  int id_ = -1;
  int size_ = -1;
  int world_size_ = -1;
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_PROCESSOR_GROUPING_HPP
