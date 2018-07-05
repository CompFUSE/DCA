// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file implements json_reader.hpp.

#include "dca/io/json/json_reader.hpp"

#include <fstream>
#include <stdexcept>

namespace dca {
namespace io {
// dca::io::

JSONReader::JSONReader()
    : current_file_name("input.json"),
      parser(),
      parse_result(parser.get_JSON_tree().result),
      my_paths(0) {}

void JSONReader::open_file(std::string file_name) {
  current_file_name = file_name;
  parse(file_name);
}

std::string JSONReader::get_path() {
  std::string path = "/";

  for (std::size_t i = 0; i < my_paths.size(); i++) {
    path = path + my_paths[i];

    if (i < my_paths.size() - 1)
      path = path + "/";
  }

  return path;
}

void JSONReader::parse(std::string& file_name_ref) {
  std::string file_name = file_name_ref;

  std::wifstream file(file_name.c_str());

  if (!file or !file.good() or file.bad()) {
    std::cout << "\n\n\tcannot open file : " << file_name << "\n";
    throw std::runtime_error(__FUNCTION__);
  }
  else {
    std::cout << "\n\n\topening file : " << file_name << "\n";
  }

  while (parser.execute(file))
    ;
}

}  // io
}  // dca
