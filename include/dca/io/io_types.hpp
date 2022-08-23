// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file defines the types of reader/writers

#ifndef DCA_IO_TYPES_HPP
#define DCA_IO_TYPES_HPP

#include <string>

namespace dca {
namespace io {
// dca::phys::

/** Enum of supported file IO types
 */
enum class IOType : int { JSON = 0, HDF5, ADIOS2 };

IOType stringToIOType(const std::string& name);

std::string toString(const IOType type);

IOType extensionToIOType(const std::string& file_name);

std::string extensionFromIOType(const IOType type);
  
}  // namespace io
}  // namespace dca

#endif  // DCA_IO_TYPES_HPP
