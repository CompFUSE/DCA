// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This class maps C++ types to MPI types.
//
// TODO: Check the MPI types.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP

#include <complex>
#include <cstdlib>
#include <mpi.h>

namespace dca {
namespace concurrency {
// dca::concurrency::

template <typename scalar_type>
class MPITypeMap {
public:
  static std::size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_PACKED;
  }
};

template <>
class MPITypeMap<bool> {
public:
  static std::size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_INT;
  }
};

template <>
class MPITypeMap<char> {
public:
  static std::size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_CHAR;
  }
};

template <>
class MPITypeMap<int> {
public:
  static std::size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_INT;
  }
};

template <>
class MPITypeMap<std::size_t> {
public:
  static std::size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_UNSIGNED_LONG;
  }
};

template <>
class MPITypeMap<float> {
public:
  static std::size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_FLOAT;
  }
};

template <>
class MPITypeMap<double> {
public:
  static std::size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_DOUBLE;
  }
};

template <>
class MPITypeMap<std::complex<float>> {
public:
  static std::size_t factor() {
    return 2;
  }

  static MPI_Datatype value() {
    return MPI_FLOAT;
  }
};

template <>
class MPITypeMap<std::complex<double>> {
public:
  static std::size_t factor() {
    return 2;
  }

  static MPI_Datatype value() {
    return MPI_DOUBLE;
  }
};

}  // concurrency
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP
