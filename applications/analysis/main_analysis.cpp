// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Main file for the DCA(+) analysis.
// Usage: ./main_analysis input_file.json

#include <iostream>
#include <string>

// Defines Concurrency, ParametersType, DcaData, and BseSolver.
#include "dca/config/analysis.hpp"

#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "comp_library/IO_library/JSON/JSON.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " input_file.json" << std::endl;
    return -1;
  }

  std::string input_file(argv[1]);

  Concurrency concurrency(argc, argv);

  // Print some info.
  if (concurrency.id() == concurrency.first()) {
    std::cout << "\nDCA(+) analysis starting.\n"
              << "MPI-world set up: " << concurrency.number_of_processors() << " processes.\n"
              << std::endl;

    dca::util::GitVersion::print();
    dca::util::Modules::print();
  }

  // Create the parameters object from the input file.
  ParametersType parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<IO::reader<IO::JSON>>(input_file);
  parameters.update_model();
  parameters.update_domains();

  // Create and initialize the DCA data object and read the output of the DCA(+) calculation.
  DcaData dca_data(parameters);
  dca_data.initialize();
  dca_data.read(parameters.get_directory() + parameters.get_output_file_name());

  // Compute the susceptibility.
  if (parameters.get_vertex_measurement_type() != NONE) {
    BseSolver analysis_obj(parameters, dca_data);
    analysis_obj.calculate_susceptibilities_2();

    if (concurrency.id() == concurrency.last()) {
      std::cout << "\nProcessor " << concurrency.id() << " is writing data " << std::endl;
      analysis_obj.write();
    }
  }

  if (concurrency.id() == concurrency.last())
    std::cout << "\nDCA(+) analysis ending.\n" << std::endl;

  return 0;
}
