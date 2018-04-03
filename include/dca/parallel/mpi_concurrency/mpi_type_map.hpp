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
// This class provides an utility for determining the size of a C++ object during MPI
// communications.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP

#include <mpi.h>

namespace dca {
namespace parallel {
// dca::parallel::

template <typename T>
struct MPITypeMap {
  static std::size_t factor() {
    return sizeof(T);
  }

  static MPI_Datatype value() {
    return MPI_CHAR;
  }
};

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP
