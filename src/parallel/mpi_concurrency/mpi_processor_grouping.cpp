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

MPIProcessorGrouping::MPIProcessorGrouping(bool (*check)()) {
  // Initialize grouping with MPI world.
  MPI_communication_ = MPI_COMM_WORLD;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_id_);

  // Check if the process has the desired qualities, or remove it from the communicator.
  const bool is_valid = check();

  MPI_Comm_split(MPI_COMM_WORLD, static_cast<int>(is_valid), 0, &MPI_communication_);

  MPI_Comm_size(MPI_communication_, &size_);

  if (is_valid) {
    MPI_Comm_rank(MPI_communication_, &id_);
    if (size_ < world_size_)
      printRemovedProcesses();
  }
  else {
    id_ = -1;
    printRemovedProcesses();
  }
}

MPIProcessorGrouping::~MPIProcessorGrouping() {
  MPI_Comm_free(&MPI_communication_);
}

bool MPIProcessorGrouping::defaultCheck() {
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

void MPIProcessorGrouping::printRemovedProcesses() const {
  // Create a communicator with removed processes plus the root of the communicator.
  MPI_Comm print_communicator;
  const int is_in_print_comm = !isValid() || id_ == 0;
  const int rank_priority = id_ == 0 ? 0 : 1;

  MPI_Comm_split(MPI_COMM_WORLD, is_in_print_comm, rank_priority, &print_communicator);
  if (!is_in_print_comm)
    return;

  int comm_id, comm_size;
  MPI_Comm_size(print_communicator, &comm_size);
  MPI_Comm_rank(print_communicator, &comm_id);

  // Get the hostname of the processors.
  char mpi_processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  const int ierr = MPI_Get_processor_name(mpi_processor_name, &name_len);
  if (ierr != MPI_SUCCESS)
    std::cerr << "MPI_Get_processor_name() failed with error " << ierr << std::endl;

  // Send hostname to root.
  std::vector<char> recvbuf(MPI_MAX_PROCESSOR_NAME * comm_size);
  MPI_Gather(mpi_processor_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, recvbuf.data(),
             MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, print_communicator);

  // Print the names.
  if (comm_id == 0) {
    std::cout
        << "\n\n**************\n"
           "\tThe following processors could not complete the required check\n"
           "\t(Default: execute a CUDA kernel) and will be removed from the communicator.\n\n";
    auto print_line = [&](const int idx) {
      std::cout << recvbuf.data() + idx * MPI_MAX_PROCESSOR_NAME << "\n";
    };

    if (!isValid())
      print_line(0);
    for (int i = 1; i < comm_size; ++i)
      print_line(i);

    std::cout << "**************\n" << std::endl;
  }

  MPI_Comm_free(&print_communicator);
}

}  // parallel
}  // dca
