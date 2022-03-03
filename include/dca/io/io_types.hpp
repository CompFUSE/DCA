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

/** Possible Flavors for G4
 *  dummy at 0 debug automatic initialization of int to 0
 *  causes a bug. That is not a good code smell.
 */
enum class IOType : int { JSON = 0, HDF5, ADIOS2, UNKNOWN };

IOType stringToIOType(const std::string& name);

std::string toString(const IOType type);

IOType extensionToIOType(const std::string& file_name);

}  // namespace io
}  // namespace dca

#endif  // DCA_PHYS_FOUR_POINT_TYPE_HPP
