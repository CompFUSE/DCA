// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef DCA_CONCURRENCY_INTERFACES_TYPE_MAP_INTERFACE_MPI_H
#define DCA_CONCURRENCY_INTERFACES_TYPE_MAP_INTERFACE_MPI_H

#include "dca/concurrency/interfaces/type_map_interface.h"
#include <mpi.h>

namespace dca {
namespace concurrency {
// dca::concurrency::

template <typename scalar_type>
class type_map_interface<MPI_LIBRARY, scalar_type> {
public:
  static size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_PACKED;
  }
};

template <>
class type_map_interface<MPI_LIBRARY, bool> {
public:
  static size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_INT;
  }
};

template <>
class type_map_interface<MPI_LIBRARY, char> {
public:
  static size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_CHAR;
  }
};

template <>
class type_map_interface<MPI_LIBRARY, int> {
public:
  static size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_INT;
  }
};

template <>
class type_map_interface<MPI_LIBRARY, size_t> {
public:
  static size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_UNSIGNED_LONG;
  }
};

template <>
class type_map_interface<MPI_LIBRARY, float> {
public:
  static size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_FLOAT;
  }
};

template <>
class type_map_interface<MPI_LIBRARY, double> {
public:
  static size_t factor() {
    return 1;
  }

  static MPI_Datatype value() {
    return MPI_DOUBLE;
  }
};

template <>
class type_map_interface<MPI_LIBRARY, std::complex<float>> {
public:
  static size_t factor() {
    return 2;
  }

  static MPI_Datatype value() {
    return MPI_FLOAT;
  }
};

template <>
class type_map_interface<MPI_LIBRARY, std::complex<double>> {
public:
  static size_t factor() {
    return 2;
  }

  static MPI_Datatype value() {
    return MPI_DOUBLE;
  }
};

}  // concurrency
}  // dca

#endif  // DCA_CONCURRENCY_INTERFACES_TYPE_MAP_INTERFACE_MPI_H
