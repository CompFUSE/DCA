// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the conversion of FourPointType to and from string.

#include "dca/io/io_types.hpp"

#include <stdexcept>

namespace dca {
namespace io {

IOType stringToIOType(const std::string& name) {
  if (name == "JSON")
    return IOType::JSON;
  else if (name == "HDF5")
    return IOType::HDF5;
  else if (name == "ADIOS2")
    return IOType::ADIOS2;
}

std::string toString(const IOType type) {
  switch (type) {
    case IOType::JSON:
      return "JSON";
    case IOType::HDF5:
      return "HDF5";
    case IOType::ADIOS2:
      return "ADIOS2";
  }
}

}  // namespace io
}  // namespace dca
