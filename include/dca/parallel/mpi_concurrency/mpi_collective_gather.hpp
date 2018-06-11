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

#include <array>
#include <map>
#include <memory>
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
  MPICollectiveGather(const std::unique_ptr<const MPIProcessorGrouping>& grouping) : grouping_(grouping) {}

  template <typename ScalarType, class Domain>
  void gather(func::function<ScalarType, Domain>& f);

  template <typename ScalarType>
  void gather(std::vector<ScalarType>& v);

private:
  template <typename ScalarType>
  void gather(ScalarType* data_ptr, int size);

  const std::unique_ptr<const MPIProcessorGrouping>& grouping_;
  std::map<int, std::array<std::vector<int>, 2>> displacements_map_;
};

template <typename ScalarType, class Domain>
void MPICollectiveGather::gather(func::function<ScalarType, Domain>& f) {
  gather(f.values(), f.size());
}

template <typename ScalarType>
void MPICollectiveGather::gather(std::vector<ScalarType>& v) {
  gather(v.data(), v.size());
}

template <typename ScalarType>
void MPICollectiveGather::gather(ScalarType* data_ptr, const int size) {
  assert(size >= 0);
  const int n_threads = grouping_->get_Nr_threads();

  std::array<std::vector<int>, 2>& entry = displacements_map_[size];
  auto& recv_count = entry[0];
  auto& displacements = entry[1];

  if (!recv_count.size()) {  // The displacements need to be computed.
    std::pair<int, int> bounds;
    recv_count.resize(n_threads);
    displacements.resize(n_threads);
    for (int id = 0; id < n_threads; ++id) {
      bounds = util::getBounds(id, n_threads, std::make_pair(0, size));
      recv_count[id] = (bounds.second - bounds.first) * MPITypeMap<ScalarType>::factor();
      displacements[id] = bounds.first * MPITypeMap<ScalarType>::factor();
    }
  }

  const auto bounds = util::getBounds(grouping_->get_id(), n_threads, std::make_pair(0, size));
  std::vector<ScalarType> sendbuff(bounds.second - bounds.first);
  std::copy_n(data_ptr + bounds.first, sendbuff.size(), sendbuff.data());

  MPI_Allgatherv(sendbuff.data(), sendbuff.size() * MPITypeMap<ScalarType>::factor(),
                 MPITypeMap<ScalarType>::value(), data_ptr, recv_count.data(), displacements.data(),
                 MPITypeMap<ScalarType>::value(), grouping_->get());
}

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_COLLECTIVE_GATHER_HPP
