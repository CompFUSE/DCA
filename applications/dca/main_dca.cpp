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

#include "dca/config/dca.hpp"
#include "dca/application/dca_loop_dispatch.hpp"
#include "dca/config/cmake_options.hpp"
// Defines Concurrency, Threading, ParametersType, DcaData, DcaLoop, and Profiler.
#include "dca/io/json/json_reader.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "dca/util/signal_handler.hpp"
#include "dca/io/writer.hpp"

int dca_main(int argc, char** argv) {
  Concurrency concurrency(argc, argv);

  try {
    // dca::util::SignalHandler<Concurrency>::init(concurrency.id() == concurrency.first());

    std::string input_file(argv[1]);

    Profiler::start();

    // Print some info.
    if (concurrency.id() == concurrency.first()) {
      dca::util::GitVersion::print();
      dca::util::Modules::print();
      dca::config::CMakeOptions::print();

#ifdef DCA_HAVE_GPU
      dca::linalg::util::printInfoDevices();
#endif  // DCA_HAVE_GPU

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

#ifdef DCA_HAVE_GPU
    dca::linalg::util::initializeMagma();
#endif  // DCA_HAVE_GPU

    // Create the parameters object from the input file.
    ParametersType parameters(dca::util::GitVersion::string(), concurrency);
    parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
    if(concurrency.id() == concurrency.first())
      std::cout << "Input read and broadcast.\n";
    parameters.update_model();
    if(concurrency.id() == concurrency.first())
      std::cout << "Model updated.\n";
    parameters.update_domains();
    if(concurrency.id() == concurrency.first())
      std::cout << "Domains updated.\n";

    dca::DistType distribution = parameters.get_g4_distribution();
    {
      switch (distribution) {
#ifdef DCA_HAVE_MPI
        case dca::DistType::BLOCKED: {
          DCALoopDispatch<dca::DistType::BLOCKED> dca_loop_dispatch;
          dca_loop_dispatch(parameters, concurrency);
        } break;
#else
        case dca::DistType::BLOCKED: {
          throw std::runtime_error(
              "Input calls for function Blocked distribution but DCA is only supports this with "
              "MPI.");
        } break;
#endif
#ifdef DCA_HAVE_MPI
        case dca::DistType::LINEAR: {
          DCALoopDispatch<dca::DistType::LINEAR> dca_loop_dispatch;
          dca_loop_dispatch(parameters, concurrency);
        } break;
#else
        case dca::DistType::LINEAR: {
          throw std::runtime_error(
              "Input calls for function Linear distribution but DCA is only supports this with "
              "MPI.");
        } break;
#endif
        case dca::DistType::NONE: {
          DCALoopDispatch<dca::DistType::NONE> dca_loop_dispatch;
          dca_loop_dispatch(parameters, concurrency);
        } break;
      }
    }
  }
  catch (const std::exception& err) {
    std::cout << "Unhandled exception in main function:\n\t" << err.what() << std::endl;
    concurrency.abort();
  }

  return 0;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " input_file.json" << std::endl;
    return -1;
  }

  return dca_main(argc, argv);
}
