// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef COMP_LIBRARY_IO_LIBRARY_CSV_CSV_WRITER_H
#define COMP_LIBRARY_IO_LIBRARY_CSV_CSV_WRITER_H

#include "comp_library/IO_library/template_writer.h"

#include <fstream>
#include <string>
#include <vector>

namespace IO {

template <>
class writer<IO::CSV> {
public:
  template <typename scalartype>
  static void execute(std::string file_name, std::vector<std::vector<scalartype>>& data);
};

template <typename scalartype>
void writer<IO::CSV>::execute(std::string file_name, std::vector<std::vector<scalartype>>& data) {
  std::ofstream myfile;
  myfile.open(file_name.c_str());

  for (size_t j = 0; j < data.size(); ++j) {
    for (size_t i = 0; i < data[j].size(); ++i) {
      myfile << data[i][j];

      if (i == data[j].size() - 1)
        myfile << "\n";
      else
        myfile << ", ";
    }
  }

  myfile.close();
}
}

#endif  // COMP_LIBRARY_IO_LIBRARY_CSV_CSV_WRITER_H
