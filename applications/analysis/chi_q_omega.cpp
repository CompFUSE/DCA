// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// Main file for generating dynamic susceptibility (Chi''_Q_omega) from a DCA++ G4
// Usage: ./chi_q_omega input_file.json

#include <iostream>
#include <string>

// Defines Concurrency, Threading, ParametersType, DcaData, and BseSolver.
#include "dca/config/analysis.hpp"
#include "dca/config/cmake_options.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/util/git_version.hpp"
#include "dca/util/modules.hpp"
#include "dca/io/writer.hpp"
#include "dca/io/reader.hpp"
#include "dca/io/io_types.hpp"

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

#ifdef DCA_HAVE_GPU
    dca::linalg::util::printInfoDevices();
#endif  // DCA_HAVE_GPU

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

#ifdef DCA_HAVE_GPU
  dca::linalg::util::initializeMagma();
#endif  // DCA_HAVE_GPU

  // Create the parameters object from the input file.
  ParametersType parameters(dca::util::GitVersion::string(), concurrency);
  parameters.read_input_and_broadcast<dca::io::JSONReader>(input_file);
  parameters.update_model();
  parameters.update_domains();

  // Create and initialize the DCA data object and read the output of the DCA(+) calculation.
  DcaDataType dca_data(parameters);
  dca_data.initialize();
  adios2::ADIOS adios;
    
#ifdef DCA_HAVE_ADIOS2
  if (dca::io::stringToIOType(parameters.get_output_format()) == dca::io::IOType::ADIOS2) {
    std::cout << "\nProcessor " << concurrency.id() << " is writing data." << std::endl;
    dca::io::Writer writer(adios, concurrency, parameters.get_output_format());
    std::string filename_bse(parameters.get_directory() + parameters.getAppropriateFilenameAnalysis());
    writer.open_file(filename_bse);

    std::string filename(parameters.get_directory() + parameters.get_filename_dca());
    dca::io::Reader reader(concurrency, parameters.get_output_format());
    reader.open_file(filename);
    auto& adios2_reader = std::get<dca::io::ADIOS2Reader<Concurrency>>(reader.getUnderlying());
    std::size_t step_count = adios2_reader.getStepCount();
    for (std::size_t i = 0; i < step_count; ++i) {
      adios2_reader.begin_step();
      writer.begin_step();
      dca_data.read(reader);
      BseSolverExt bse_solver_ext(parameters, dca_data);
      bse_solver_ext.calculateSusceptibilities();
      bse_solver_ext.write(writer);
      adios2_reader.end_step();
      writer.end_step();
    }
  }
#endif
  std::cout << "\nFinish time: " << dca::util::print_time() << "\n" << std::endl;
return 0;
}
