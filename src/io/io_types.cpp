// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This file implements the conversion of IOType to and from string.

#include "dca/io/io_types.hpp"
#include <iostream>

namespace dca {
namespace io {

IOType stringToIOType(const std::string& name) {
  if (name == "JSON")
    return IOType::JSON;
  else if (name == "HDF5")
    return IOType::HDF5;
  else if (name == "ADIOS2")
    return IOType::ADIOS2;
  else
    throw std::runtime_error(name + " is an unknown file format identifier!");
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

IOType extensionToIOType(const std::string& file_name) {
  // look for the first extension
  std::size_t extension_start = file_name.rfind('.');
  std::string extension{file_name.substr(extension_start + 1, file_name.size())};
  if (extension == "bp")
    return IOType::ADIOS2;
  else if (extension == "hdf" || extension == "hdf5" || extension == "tmp")
    return IOType::HDF5;
  else if (extension == "json")
    return IOType::JSON;
  else
    throw std::runtime_error("extension " + extension + " identifies an unknown file format\n"
			     + "known extensions (bp|hdf|hdf5|tmp|json)");
}

std::string extensionFromIOType(const IOType type) {
  switch (type) {
    case IOType::JSON:
      return ".json";
    case IOType::HDF5:
      return ".hdf5";
    case IOType::ADIOS2:
      return ".bp";
  }
}

  
}  // namespace io
}  // namespace dca
