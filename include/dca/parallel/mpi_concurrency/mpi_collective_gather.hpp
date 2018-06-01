// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class gathers data computed on different processes stored in a portion of a contiguous
// buffer determined by dca::parallel::util::getBounds.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_GATHER_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_GATHER_HPP

#include <vector>

#include <mpi.h>

#include "dca/function/function.hpp"
#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"
#include "dca/parallel/mpi_concurrency/mpi_type_map.hpp"
#include "dca/parallel/util/get_bounds.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

class MPICollectiveGather {
public:
  MPICollectiveGather(const MPIProcessorGrouping& grouping) : grouping_(grouping) {}

  template <typename ScalarType, class Domain>
  void gather(func::function<ScalarType, Domain>& f) const;

  template <typename ScalarType>
  void gather(std::vector<ScalarType>& v) const;

private:
  template <typename ScalarType>
  void gather(ScalarType* data_ptr, int size) const;

  const MPIProcessorGrouping& grouping_;
};

template <typename ScalarType, class Domain>
void MPICollectiveGather::gather(func::function<ScalarType, Domain>& f) const {
  gather(f.values(), f.size());
}

template <typename ScalarType>
void MPICollectiveGather::gather(std::vector<ScalarType>& v) const {
  gather(v.data(), v.size());
}

template <typename ScalarType>
void MPICollectiveGather::gather(ScalarType* data_ptr, const int size) const {
  assert(size >= 0);

  const int n_threads = grouping_.get_Nr_threads();
  std::vector<int> displs(n_threads);
  std::vector<int> recvcount(n_threads);

  std::pair<int, int> bounds;
  for (int id = 0; id < n_threads; ++id) {
    bounds = util::getBounds(id, n_threads, std::make_pair(0, size));
    displs[id] = bounds.first * MPITypeMap<ScalarType>::factor();
    recvcount[id] = (bounds.second - bounds.first) * MPITypeMap<ScalarType>::factor();
  }

  bounds = util::getBounds(grouping_.get_id(), n_threads, std::make_pair(0, size));
  std::vector<ScalarType> sendbuff(bounds.second - bounds.first);
  std::copy_n(data_ptr + bounds.first, sendbuff.size(), sendbuff.data());

  MPI_Allgatherv(sendbuff.data(), sendbuff.size() * MPITypeMap<ScalarType>::factor(),
                 MPITypeMap<ScalarType>::value(), data_ptr, recvcount.data(), displs.data(),
                 MPITypeMap<ScalarType>::value(), grouping_.get());
}

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_GATHER_HPP
