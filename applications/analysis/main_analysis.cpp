// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// Main file for the DCA(+) analysis.
// Usage: ./main_analysis input_file.json

#include <iostream>
#include <string>

// Defines Concurrency, Threading, ParametersType, DcaData, and BseSolver.
#include "dca/config/analysis.hpp"
#include "dca/config/cmake_options.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " input_file.json" << std::endl;
    return -1;
  }

  std::string input_file(argv[1]);

  Concurrency concurrency(argc, argv);

  // Print some info.
  if (concurrency.id() == concurrency.first()) {
    dca::util::GitVersion::print();
    dca::util::Modules::print();
    dca::config::CMakeOptions::print();

    std::cout
        << "\n"
        << "********************************************************************************\n"
        << "**********                      DCA(+) Analysis                       **********\n"
        << "********************************************************************************\n"
        << "\n"
        << "Start time : " << dca::util::print_time() << "\n"
        << "\n"
        << "MPI-world set up: " << concurrency.number_of_processors() << " processes."
        << "\n"
        << std::endl;
  }

  // Create the parameters object from the input file.
  ParametersType parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
  parameters.update_model();
  parameters.update_domains();

  // Create and initialize the DCA data object and read the output of the DCA(+) calculation.
  DcaDataType dca_data(parameters);
  dca_data.initialize();
  dca_data.read(parameters.get_directory() + parameters.get_filename_dca());

  BseSolverType bse_solver(parameters, dca_data);
  bse_solver.calculateSusceptibilities();

  if (concurrency.id() == concurrency.first()) {
    std::cout << "\nProcessor " << concurrency.id() << " is writing data." << std::endl;
    bse_solver.write();

    std::cout << "\nFinish time: " << dca::util::print_time() << "\n" << std::endl;
  }

  return 0;
}
