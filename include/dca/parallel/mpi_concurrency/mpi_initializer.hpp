// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class initializes the standard Message Passing Interface (MPI).

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_INITIALIZER
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_INITIALIZER

namespace dca {
namespace parallel {
// dca::parallel::

class MPIInitializer {
protected:
  MPIInitializer(int argc, char** argv);

  ~MPIInitializer();

public:
  // Aborts all processes.
  void abort() const;
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_INITIALIZER
