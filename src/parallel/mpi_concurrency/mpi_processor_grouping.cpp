// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements mpi_processor_grouping.hpp.

#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"

#include <vector>

#include "dca/linalg/device_type.hpp"
#include "dca/linalg/util/memory.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

MPIProcessorGrouping::MPIProcessorGrouping() {
  // Initialize grouping with MPI world.
  MPI_communication_ = MPI_COMM_WORLD;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_threads_);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_id_);

  // Check if the process has the desired qualities, or remove it from the communicator.
  const bool is_valid_ = testValidity();

  std::vector<int> validity_input(nr_threads_, false);
  std::vector<int> validity_output(nr_threads_, false);
  validity_input[world_id_] = is_valid_;
  MPI_Allreduce(validity_input.data(), validity_output.data(), validity_input.size(), MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);

  std::vector<int> valid_ids;
  for (int id = 0; id < validity_output.size(); ++id)
    if (validity_output[id])
      valid_ids.push_back(id);
  if (!valid_ids.size())
    throw(std::logic_error("MPI grouping is empty.\n"));

  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group_incl(world_group, valid_ids.size(), valid_ids.data(), &MPI_group_);
  MPI_Comm_create_group(MPI_COMM_WORLD, MPI_group_, 1, &MPI_communication_);
  MPI_Group_free(&world_group);

  if (is_valid_) {
    MPI_Comm_size(MPI_communication_, &nr_threads_);
    MPI_Comm_rank(MPI_communication_, &id_);
  }
  else {
    id_ = nr_threads_ = -1;
  }
}

MPIProcessorGrouping::~MPIProcessorGrouping() {
  MPI_Comm_free(&MPI_communication_);
  MPI_Group_free(&MPI_group_);
}

bool MPIProcessorGrouping::testValidity() const {
#ifdef DCA_HAVE_CUDA
  constexpr linalg::DeviceType device = linalg::GPU;
#else
  constexpr linalg::DeviceType device = linalg::CPU;
#endif  // DCA_HAVE_CUDA

  try {
    double* ptr = nullptr;
    linalg::util::Memory<device>::allocate(ptr, 32);
    linalg::util::Memory<device>::deallocate(ptr);

    return true;
  }
  catch (...) {
    return false;
  }
}

}  // parallel
}  // dca
