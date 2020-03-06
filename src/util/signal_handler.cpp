// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This singleton class organizes the handling of signals.

#include <csignal>
#include <iostream>

#include "dca/util/signal_handler.hpp"

namespace dca {
namespace util {
// dca::util::

void SignalHandler::init(bool verbose) {
  verbose_ = verbose;

  signal(SIGINT, sigintHandler);
}

void SignalHandler::sigintHandler(int signum) {
  if (verbose_)
    std::cerr << "Interrupt signal (" << signum << ") received.\n";

  for (auto file : file_ptrs_)
    file->close_file();

  exit(signum);
}

void SignalHandler::registerFile(io::HDF5Writer& writer) {
  file_ptrs_.push_back(&writer);
}

}  // namespace util
}  // namespace dca
