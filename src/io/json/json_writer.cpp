// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements json_writer.hpp.

#include "dca/io/json/json_writer.hpp"
#include <fstream>

namespace dca {
namespace io {
// dca::io::

JSONWriter::JSONWriter() : file_name(""), path(""), elements_in_group(0) {
  ss << std::fixed;
  ss.precision(16);
}

std::stringstream& JSONWriter::open_file(const std::string& my_file_name, const bool /*overwrite*/) {
  file_name = my_file_name;

  ss << "{\n";

  elements_in_group.push_back(0);

  return ss;
}

void JSONWriter::close_file() {
  ss << "\n}";

  std::ofstream of(&file_name[0]);

  of << ss.str();

  of.flush();
  of.close();
}

void JSONWriter::open_group(const std::string& name) {
  if (elements_in_group.back() != 0)
    ss << ",\n\n";

  ss << get_path() << "\"" << name << "\""
     << " : \n";

  ss << get_path() << "{\n";

  elements_in_group.push_back(0);
}

void JSONWriter::close_group() {
  elements_in_group.pop_back();

  ss << "\n" << get_path() << "}";

  elements_in_group.back() += 1;
}

std::string JSONWriter::get_path() {
  std::stringstream ss;
  for (size_t i = 0; i < elements_in_group.size(); i++)
    ss << "\t";

  return ss.str();
}

void JSONWriter::execute(const std::string& name, const std::string& value) {
  if (elements_in_group.back() != 0)
    ss << ",\n";

  ss << get_path() << "\"" << name << "\" : \"" << value << "\"";

  elements_in_group.back() += 1;
}

void JSONWriter::execute(const std::string& name, const std::vector<std::string>& value) {
  if (elements_in_group.back() != 0)
    ss << ",\n";

  ss << get_path() << "\"" << name << "\" : [";

  for (size_t i = 0; i < value.size(); i++) {
    ss << "\"" << value[i] << "\"";

    if (i == value.size() - 1)
      ss << "]";
    else
      ss << ", ";
  }

  elements_in_group.back() += 1;
}

}  // io
}  // dca
