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
struct MPITypeMap {
  template <typename = std::enable_if_t<std::is_enum<T>::value>>
  static MPI_Datatype value() {
    return MPITypeMap<std::underlying_type_t<T>>::value();
  }
};

template <>
struct MPITypeMap<bool> {
  static MPI_Datatype value() {
    return MPI_CXX_BOOL;
  }
};

template <>
struct MPITypeMap<char> {
  static MPI_Datatype value() {
    return MPI_CHAR;
  }
};

template <>
struct MPITypeMap<std::uint8_t> {
  static MPI_Datatype value() {
    return MPI_UNSIGNED_CHAR;
  }
};


template <>
struct MPITypeMap<int> {
  static MPI_Datatype value() {
    return MPI_INT;
  }
};

template <>
struct MPITypeMap<std::size_t> {
  static MPI_Datatype value() {
    return MPI_UNSIGNED_LONG;
  }
};

template <>
struct MPITypeMap<long long int> {
  static MPI_Datatype value() {
    return MPI_LONG_LONG_INT;
  }
};

template <>
struct MPITypeMap<float> {
  static MPI_Datatype value() {
    return MPI_FLOAT;
  }
};

template <>
struct MPITypeMap<double> {
  static MPI_Datatype value() {
    return MPI_DOUBLE;
  }
};

template <>
struct MPITypeMap<std::complex<float>> {
  static MPI_Datatype value() {
    return MPI_COMPLEX;
  }
};

template <>
struct MPITypeMap<std::complex<double>> {
  static MPI_Datatype value() {
    return MPI_DOUBLE_COMPLEX;
  }
};

}  // namespace parallel
}  // namespace dca

#endif  // DCA_PARALLEL_MPI_CONCURRENCY_MPI_TYPE_MAP_HPP
