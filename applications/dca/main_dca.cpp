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
// Main file for the DCA(+) calculation.
// Usage: ./main_dca input_file.json

#include <string>
#include <iostream>

#include "dca/config/cmake_options.hpp"
// Defines Concurrency, Threading, ParametersType, DcaData, DcaLoop, and Profiler.
#include "dca/config/dca.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " input_file.json" << std::endl;
    return -1;
  }

  Concurrency concurrency(argc, argv);

  try {
    std::string input_file(argv[1]);

    Profiler::start();

    // Print some info.
    if (concurrency.id() == concurrency.first()) {
      dca::util::GitVersion::print();
      dca::util::Modules::print();
      dca::config::CMakeOptions::print();

#ifdef DCA_WITH_CUDA
      dca::linalg::util::printInfoDevices();
#endif  // DCA_WITH_CUDA

      std::cout
          << "\n"
          << "********************************************************************************\n"
          << "**********                     DCA(+) Calculation                     **********\n"
          << "********************************************************************************\n"
          << "\n"
          << "Start time : " << dca::util::print_time() << "\n"
          << "\n"
          << "MPI-world set up: " << concurrency.number_of_processors() << " processes."
          << "\n"
          << std::endl;
    }

#ifdef DCA_WITH_CUDA
    dca::linalg::util::initializeMagma();
#endif  // DCA_WITH_CUDA

    // Create the parameters object from the input file.
    ParametersType parameters(dca::util::GitVersion::string(), concurrency);
    parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
    parameters.update_model();
    parameters.update_domains();

    // Create and initialize the DCA data object.
    DcaDataType dca_data(parameters);
    dca_data.initialize();

    DcaLoopType dca_loop(parameters, dca_data, concurrency);

    {
      Profiler profiler(__FUNCTION__, __FILE__, __LINE__);

      dca_loop.initialize();
      dca_loop.execute();
      dca_loop.finalize();
    }

    Profiler::stop(concurrency, parameters.get_filename_profiling());

    if (concurrency.id() == concurrency.first()) {
      std::cout << "\nProcessor " << concurrency.id() << " is writing data." << std::endl;
      dca_loop.write();

      std::cout << "\nFinish time: " << dca::util::print_time() << "\n" << std::endl;
    }
  }
  catch (const std::exception& err) {
    std::cout << "Unhandled exception in main function:\n\t" << err.what();
    concurrency.abort();
  }

  return 0;
}
