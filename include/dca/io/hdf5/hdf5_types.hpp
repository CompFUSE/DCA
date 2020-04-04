// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// HDF5 types.

#ifndef DCA_IO_HDF5_HDF5_TYPES_HPP
#define DCA_IO_HDF5_HDF5_TYPES_HPP

#include <complex>
#include "H5Cpp.h"

namespace dca {
namespace io {
// dca::io::

template <typename HDF5_type>
class HDF5_TYPE {};

template <>
class HDF5_TYPE<bool> {
public:
  static hid_t get() {
    return H5T_NATIVE_HBOOL;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_HBOOL;
  }
};

template <>
class HDF5_TYPE<int> {
public:
  static hid_t get() {
    return H5T_NATIVE_INT;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_INT;
  }
};

template <>
class HDF5_TYPE<unsigned int> {
public:
  static hid_t get() {
    return H5T_NATIVE_UINT;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_UINT;
  }
};

template <>
class HDF5_TYPE<long> {
public:
  static hid_t get() {
    return H5T_NATIVE_LONG;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_LONG;
  }
};

template <>
class HDF5_TYPE<unsigned long> {
public:
  static hid_t get() {
    return H5T_NATIVE_ULONG;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_ULONG;
  }
};

template <>
class HDF5_TYPE<long long> {
public:
  static hid_t get() {
    return H5T_NATIVE_LLONG;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_LLONG;
  }
};

template <>
class HDF5_TYPE<unsigned long long> {
public:
  static hid_t get() {
    return H5T_NATIVE_ULLONG;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_ULLONG;
  }
};

template <>
class HDF5_TYPE<char> {
public:
  static hid_t get() {
    return H5T_NATIVE_CHAR;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_CHAR;
  }
};

template <>
class HDF5_TYPE<float> {
public:
  static hid_t get() {
    return H5T_NATIVE_FLOAT;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_FLOAT;
  }
};

template <>
class HDF5_TYPE<double> {
public:
  static hid_t get() {
    return H5T_NATIVE_DOUBLE;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_DOUBLE;
  }
};

template <>
class HDF5_TYPE<unsigned char> {
public:
  static hid_t get() {
    return H5T_NATIVE_UCHAR;
  }

  static H5::PredType get_PredType() {
    return H5::PredType::NATIVE_UCHAR;
  }
};

template <class T>
class HDF5_TYPE<std::complex<T>> {
public:
  static H5::CompType get_PredType() {
    auto build_type = []() {
      H5::CompType complex_data_type(2 * sizeof(T));
      complex_data_type.insertMember("r", 0, HDF5_TYPE<T>::get_PredType());
      complex_data_type.insertMember("i", sizeof(T), HDF5_TYPE<T>::get_PredType());
      return complex_data_type;
    };
    static H5::CompType type = build_type();
    return type;
  }
};

}  // namespace io
}  // namespace dca

#endif  // DCA_IO_HDF5_HDF5_TYPES_HPP
