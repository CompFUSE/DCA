// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Norbert podhorszki (pnorbert@ornl.gov)
//
// This file implements adios2_reader.hpp.

#include "dca/io/adios2/adios2_reader.hpp"

#include <fstream>
#include <stdexcept>

namespace dca {
namespace io {
// dca::io::

ADIOS2Reader::ADIOS2Reader(const std::string& config, bool verbose)
    : adios_(adios2::ADIOS(config)), verbose_(verbose) {}

ADIOS2Reader::~ADIOS2Reader() {
  if (file_)
    close_file();
}

void ADIOS2Reader::open_file(const std::string& file_name) {
  if (verbose_) {
    std::cout << "\n\n\topening ADIOS2 file : " << file_name << "\n";
  }

  io_name_ = file_name;
  file_name_ = file_name;
  io_ = adios_.DeclareIO(io_name_);
  file_ = io_.Open(file_name_, adios2::Mode::Read);
}

void ADIOS2Reader::close_file() {
  if (file_) {
    file_.Close();
    adios_.RemoveIO(io_name_);
  }
}

void ADIOS2Reader::open_group(const std::string& name) {
  my_paths_.push_back(name);
}

void ADIOS2Reader::close_group() {
  my_paths_.pop_back();
}

std::string ADIOS2Reader::get_path(const std::string& name) {
  std::string path = "/";

  for (size_t i = 0; i < my_paths_.size(); i++) {
    path += my_paths_[i];

    if (i < my_paths_.size() - 1)
      path += "/";
  }

  if (!name.empty()) {
    path += name;
  }

  return path;
}

bool ADIOS2Reader::execute(const std::string& name, std::string& value) {
  std::string full_name = get_path(name);
  if (!exists(full_name)) {
    return false;
  }

  /* ADIOS will resize string to match size of incoming data */
  file_.Get<std::string>(full_name, value, adios2::Mode::Sync);

  return true;
}

bool ADIOS2Reader::execute(const std::string& name, std::vector<std::string>& value) {
  std::string full_name = get_path(name);
  if (!exists(full_name)) {
    return false;
  }

  value.resize(1);

  /* ADIOS will resize string to match size of incoming data */
  file_.Get<std::string>(full_name, value[0], adios2::Mode::Sync);

  return true;
}

bool ADIOS2Reader::exists(const std::string& name) const {
  std::string varType = io_.VariableType(name);
  return !varType.empty();
}

}  // namespace io
}  // namespace dca
