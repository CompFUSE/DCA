// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file implements timer.hpp.

#include "dca/util/timer.hpp"
#include "dca/util/print_time.hpp"
#include <iostream>

namespace dca {
namespace util {
// dca::util::

Timer::Timer(const std::string& name, const bool output)
    : name_(name), output_(output), start_(std::chrono::system_clock::now()) {
  if (output_)
    std::cout << name_ << " start:    " << print_time(start_) << std::endl;
}

Timer::~Timer() {
  if (output_) {
    const auto end(std::chrono::system_clock::now());
    const std::chrono::duration<double> elapsed_seconds(end - start_);

    std::cout << name_ << " end:      " << print_time(end) << "\n"
              << name_ << " duration: " << elapsed_seconds.count() << " s"
              << "\n"
              << std::endl;
  }
}

}  // util
}  // dca
