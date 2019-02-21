// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides a utility function that executes a system call and returns the standard
// output.

#ifndef DCA_UTIL_GET_STDOUT_FROM_COMMAND_HPP
#define DCA_UTIL_GET_STDOUT_FROM_COMMAND_HPP

#include <string>

namespace dca {
namespace util {

std::string getStdoutFromCommand(std::string cmd);

}  // namespace util
}  // namespace dca

#endif  // DCA_UTIL_GET_STDOUT_FROM_COMMAND_HPP
