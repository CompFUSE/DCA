// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
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

  signal(SIGABRT, handle);
  signal(SIGINT, handle);
  signal(SIGFPE, handle);
  signal(SIGILL, handle);
  signal(SIGSEGV, handle);
  signal(SIGTERM, handle);
  signal(SIGUSR2, handle);  // Summit out of time signal.
}

void SignalHandler::handle(int signum) {
  if (verbose_)
    std::cerr << "Received signal (" << signum << ") received." << std::endl;

  for (auto file_ptr : file_ptrs_) {
    auto file = file_ptr.lock();
    if (file)
      file->close_file();
  }

  exit(signum);
}

void SignalHandler::registerFile(const std::shared_ptr<io::Writer>& writer) {
  file_ptrs_.emplace_back(writer);
}

}  // namespace util
}  // namespace dca
