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
#include "dca/config/mc_options.hpp"
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
        << "**********                      DCA(+) S(q,omega)  Calculation        **********\n"
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

#ifdef DCA_HAVE_ADIOS2
  adios2::ADIOS adios;

  if (dca::io::stringToIOType(parameters.get_output_format()) == dca::io::IOType::ADIOS2) {
    int rank = concurrency.id();
    std::cout << "\nProcessor " << concurrency.id() << " is writing data." << std::endl;
    dca::io::Writer writer(concurrency, parameters.get_output_format(), true);
    std::string filename_bse(parameters.get_directory() + parameters.getAppropriateFilenameAnalysis());
    writer.open_file(filename_bse);

    std::string filename(parameters.get_directory() + parameters.get_filename_dca());
    dca::io::Reader reader(concurrency, parameters.get_output_format());
    reader.open_file(filename);
    auto& adios2_reader = std::get<dca::io::ADIOS2Reader<Concurrency>>(reader.getUnderlying());
    std::size_t step_count = adios2_reader.getStepCount();

    int ranks = concurrency.get_size();

    int current_step = 0;

    while (true) {
      if (current_step >= step_count - 1)
        break;      
      writer.begin_step();
      int got_step_num = -1;
      for (int ir = 0; ir < ranks; ++ir) {
	// cut off the last weird step dca is saving
        if (current_step >= step_count - 1)
          goto end_steps;
        adios2_reader.begin_step();
        if (rank == ir) {
          dca_data.read(reader);
          got_step_num = current_step;
        }
        adios2_reader.end_step();
        ++current_step;
      }
      if (got_step_num >= 0) {
        std::cout << "Rank " << rank << " working on step " << got_step_num << '\n';
        BseSolverExt bse_solver_ext(parameters, dca_data);
        bse_solver_ext.calculateSusceptibilities();
        writer.open_group("step_" + std::to_string(got_step_num));
        bse_solver_ext.write(writer);
        writer.close_group();
        got_step_num = -1;
      }
    end_steps:
      writer.end_step();
    }
    writer.begin_step();
    if (concurrency.id() == concurrency.first())
      parameters.write(writer);
    writer.end_step();
    writer.close_file();
  }
#endif
  std::cout << "\nFinish time: " << dca::util::print_time() << "\n" << std::endl;
  return 0;
}
