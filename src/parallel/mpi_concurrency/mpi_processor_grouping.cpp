// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements mpi_processor_grouping.hpp.

#include "dca/parallel/mpi_concurrency/mpi_processor_grouping.hpp"

#include <iostream>
#include <unistd.h>

#include "dca/parallel/mpi_concurrency/kernel_test.hpp"

namespace dca {
namespace parallel {
// dca::parallel::

MPIProcessorGrouping::MPIProcessorGrouping() : MPIProcessorGrouping(defaultTest) {}

MPIProcessorGrouping::MPIProcessorGrouping(bool (*test)()) {
  // Initialize grouping with MPI world.
  MPI_communication_ = MPI_COMM_WORLD;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_id_);

  // Check if the process has the desired qualities, or remove it from the communicator.
  const bool is_valid = test();

  std::vector<char> validity_input(world_size_, 0);
  std::vector<char> validity_output(world_size_, 0);
  validity_input[world_id_] = is_valid;
  MPI_Allreduce(validity_input.data(), validity_output.data(), validity_input.size(), MPI_CHAR,
                MPI_MAX, MPI_COMM_WORLD);

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

  if (is_valid) {
    MPI_Comm_size(MPI_communication_, &size_);
    MPI_Comm_rank(MPI_communication_, &id_);
  }
  else {
    id_ = size_ = -1;
  }

  if (world_size_ > size_)
    printRemovedProcesses(valid_ids);
}

MPIProcessorGrouping::~MPIProcessorGrouping() {
  if (isValid()) {
    MPI_Comm_free(&MPI_communication_);
    MPI_Group_free(&MPI_group_);
  }
}

bool MPIProcessorGrouping::defaultTest() {
#ifdef DCA_HAVE_CUDA
  try {
    return kernelTest();
  }
  catch (...) {
    return false;
  }
#else
  return true;
#endif  // DCA_HAVE_CUDA
}

void MPIProcessorGrouping::printRemovedProcesses(const std::vector<int>& valid_ids) const {
  constexpr std::size_t len = 50;

  // Create a communicator with removed processes plus the root of the communicator.
  std::vector<int> removed_ids;
  if (valid_ids.size())
    removed_ids.push_back(valid_ids[0]);
  for (int id = 0, confornt_idx = 0; id < world_size_; ++id) {
    if (id != valid_ids[confornt_idx])
      removed_ids.push_back(id);
    else
      ++confornt_idx;
  }

  if (isValid() && removed_ids[0] != world_id_)
    return;

  MPI_Group world_group, removed_group;
  MPI_Comm removed_comm;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group_incl(world_group, removed_ids.size(), removed_ids.data(), &removed_group);
  MPI_Comm_create_group(MPI_COMM_WORLD, removed_group, 2, &removed_comm);
  MPI_Group_free(&world_group);

  // Send hostname to root.
  std::vector<char> recv_buffer;
  if (world_id_ == removed_ids[0])
    recv_buffer.resize(len * removed_ids.size());

  char hostname[len];
  gethostname(hostname, len);
  MPI_Gather(hostname, len, MPI_CHAR, recv_buffer.data(), len, MPI_CHAR, 0, removed_comm);

  // Print the information.
  if (world_id_ == removed_ids[0]) {
    std::cout << " \n\n\n********* Invalid processes location *********\n";
    auto print_line = [&](const int idx) { std::cout << recv_buffer.data() + idx * len << "\n"; };
    if (!isValid())
      print_line(0);
    for (int i = 1; i < removed_ids.size(); ++i)
      print_line(i);

    std::cout << "**********************************************\n" << std::endl;
  }

  MPI_Comm_free(&removed_comm);
  MPI_Group_free(&removed_group);
}

}  // parallel
}  // dca
