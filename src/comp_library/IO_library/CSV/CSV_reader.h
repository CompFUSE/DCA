// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//
// Description

#ifndef COMP_LIBRARY_IO_LIBRARY_CSV_CSV_READER_H
#define COMP_LIBRARY_IO_LIBRARY_CSV_CSV_READER_H

#include "comp_library/IO_library/template_reader.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace IO {

template <>
class reader<IO::CSV> {
public:
  template <typename scalartype>
  static void execute(std::string file_name, std::vector<std::vector<scalartype>>& data);
};

template <typename scalartype>
void reader<IO::CSV>::execute(std::string file_name, std::vector<std::vector<scalartype>>& data) {
  std::filebuf fb;

  if (fb.open(file_name.c_str(), std::ios::in)) {
    std::istream myfile(&fb);

    std::string row;

    data.resize(0);
    while (std::getline(myfile, row)) {
      data.push_back(std::vector<scalartype>());

      std::istringstream tokenS(row);
      std::string token;

      while (std::getline(tokenS, token, ',')) {
        std::istringstream valueS(token);

        valueS.imbue(myfile.getloc());

        scalartype value;

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
}

#endif  // COMP_LIBRARY_IO_LIBRARY_CSV_CSV_READER_H
