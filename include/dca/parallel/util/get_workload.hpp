// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a function that computes the local work size given local id and total amount
// of work.

#ifndef DCA_PARALLEL_UTIL_GET_WORKLOAD_HPP
#define DCA_PARALLEL_UTIL_GET_WORKLOAD_HPP

#include <cassert>

namespace dca {
namespace parallel {
namespace util {
// dca::parallel::util::

inline int getWorkload(const unsigned int total_work, const unsigned int n_workers,
                       const unsigned int id) {
  int work = total_work / n_workers;
  if (id < (total_work % n_workers))
    ++work;
  return work;
}

template <class Concurrency>
int getWorkload(const unsigned int total_work, const Concurrency& concurrency) {
  return getWorkload(total_work, concurrency.number_of_processors(), concurrency.id());
}

template <class Concurrency>
int getWorkload(const unsigned int total_work, const unsigned int n_local_workers,
                const unsigned int local_id, const Concurrency& concurrency) {
  assert(local_id < n_local_workers);
  const int local_work = getWorkload(total_work, concurrency);

  return getWorkload(local_work, n_local_workers, local_id);
}

inline void getComputeRange(const int& my_rank, const int& mpi_size,
                            const uint64_t& total_G4_size, uint64_t& start, uint64_t& end) {
  uint64_t offset = 0;
  // check if originally flattened one-dimensional G4 array can be equally (up to 0) distributed across ranks
  // if balanced, each rank has same amount of elements to compute
  // if not, ranks with (rank_id < nb_more_work_ranks) has to compute 1 more element than other ranks
  bool balanced = (total_G4_size % static_cast<uint64_t>(mpi_size) == 0);
  uint64_t local_work = total_G4_size / static_cast<uint64_t>(mpi_size);

  if(balanced) {
        offset = static_cast<uint64_t>(my_rank)  * local_work;
        end  = offset + local_work;
  }
  else {
    int nb_more_work_ranks = total_G4_size % static_cast<uint64_t>(mpi_size);

    if (my_rank < nb_more_work_ranks) {
      offset = static_cast<uint64_t>(my_rank) * (local_work + 1);
      end = offset + (local_work + 1);
    } else {
        offset = nb_more_work_ranks * (local_work + 1) +
                 (static_cast<uint64_t>(my_rank) - nb_more_work_ranks) * local_work;
        end = offset + local_work;
      }
    }
    start = offset;
}

}  // util
}  // parallel
}  // dca

#endif  // DCA_PARALLEL_UTIL_GET_WORKLOAD_HPP
