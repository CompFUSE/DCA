// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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
