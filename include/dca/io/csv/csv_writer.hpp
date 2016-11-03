// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// CSV writer.

#ifndef DCA_IO_CSV_CSV_WRITER_HPP
#define DCA_IO_CSV_CSV_WRITER_HPP

#include <fstream>
#include <string>
#include <vector>

namespace dca {
namespace io {
// dca::io::

class CSVWriter {
public:
  template <typename ScalarType>
  static void execute(const std::string& file_name, const std::vector<std::vector<ScalarType>>& data);
};

template <typename ScalarType>
void CSVWriter::execute(const std::string& file_name,
                        const std::vector<std::vector<ScalarType>>& data) {
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

}  // io
}  // dca

#endif  // DCA_IO_CSV_CSV_WRITER_HPP
