// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// JSON reader.

#include "dca/io/json/json_reader.hpp"

#include <fstream>

namespace dca::io {

JSONReader::JSONReader(bool verbose) : verbose_(verbose) {
  open_groups_.push(&root_);
}

void JSONReader::open_file(const std::string& filename) {
  std::ifstream inp(filename);

  if (!inp) {
    throw(std::runtime_error("Can not open file " + filename));
  }

  std::stringstream stream;
  stream << inp.rdbuf();
  inp.close();

  try {
    root_.read(stream);
  }
  catch (const std::logic_error& err) {
    throw(std::logic_error("File " + filename + ":\n" + err.what()));
  }
}

void JSONReader::close_file() noexcept {
  root_.clear();
  while (!open_groups_.empty())
    open_groups_.pop();
}

bool JSONReader::open_group(const std::string& name) noexcept {
  details::JSONGroup* new_group = nullptr;
  if (open_groups_.top())
    new_group = open_groups_.top()->getGroup(name);

  // TODO maybe: process error here.
  //  if (!new_group)
  //    throw(std::logic_error("Group " + name + " does not exist"));
  open_groups_.push(new_group);
  return static_cast<bool>(new_group);
}

bool JSONReader::close_group() noexcept {
  if (open_groups_.size() == 1) {
    //    throw(std::logic_error("Can't close root group."));
    return false;
  }

  open_groups_.pop();
  return true;
}

}  // namespace dca::io
