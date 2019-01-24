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
// This class maps C++ types to MPI types.
//
// TODO: Check the MPI types.

#ifndef DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP
#define DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP

#include <complex>
#include <cstdlib>
#include <type_traits>
#include <mpi.h>

namespace dca {
namespace parallel {
// dca::parallel::

template <typename T>
class MPITypeMap {
public:
  static constexpr std::size_t factor() {
    return 1;
  }

  template <typename = std::enable_if_t<std::is_enum<T>::value>>
  static MPI_Datatype value() {
    return MPITypeMap<std::underlying_type_t<T>>::value();
  }
};

template <>
class MPITypeMap<bool> {
public:
  static std::size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_CXX_BOOL;
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
class MPITypeMap<long long int> {
public:
  static std::size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_LONG_LONG_INT;
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

}  // parallel
}  // dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP
