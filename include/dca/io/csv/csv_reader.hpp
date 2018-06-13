// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// CSV reader.

#ifndef DCA_IO_CSV_CSV_READER_HPP
#define DCA_IO_CSV_CSV_READER_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace dca {
namespace io {
// dca::io::

class CSVReader {
public:
  template <typename ScalarType>
  static void execute(const std::string& file_name, std::vector<std::vector<ScalarType>>& data);
};

template <typename ScalarType>
void CSVReader::execute(const std::string& file_name, std::vector<std::vector<ScalarType>>& data) {
  std::filebuf fb;

  if (fb.open(file_name.c_str(), std::ios::in)) {
    std::istream myfile(&fb);

    std::string row;

    data.resize(0);
    while (std::getline(myfile, row)) {
      data.push_back(std::vector<ScalarType>());

      std::istringstream tokenS(row);
      std::string token;

      while (std::getline(tokenS, token, ',')) {
        std::istringstream valueS(token);

        valueS.imbue(myfile.getloc());

        ScalarType value;

        if (valueS >> value)
          data.back().push_back(value);
      }
    }

    fb.close();
  }
  else {
    std::cout << "\n\n\t " << file_name << " can not be read !!! \n\n\t";
    throw std::logic_error(__FUNCTION__);
  }
}

}  // io
}  // dca

#endif  // DCA_IO_CSV_CSV_READER_HPP
