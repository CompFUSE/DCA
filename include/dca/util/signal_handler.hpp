// Copyright (C) 2019 ETH Zurich
// Copyright (C) 2019 UT-Battelle, LLC
// All rights reserved.
// See LICENSE.txt for terms of usage./
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This singleton class organizes the handling of signals.

#ifndef DCA_UTIL_SIGNAL_HANDLER_HPP
#define DCA_UTIL_SIGNAL_HANDLER_HPP

#include <memory>
#include <vector>

#include "dca/io/writer.hpp"

namespace dca {
namespace util {
// dca::util::

template<class Concurrency>
class SignalHandler{
public:
    static void init(bool verbose = false);

  static void registerFile(const std::shared_ptr<io::Writer<Concurrency>>& writer);

private:
    static void handle(int signum);

    static inline bool verbose_;
  static inline std::vector<std::weak_ptr<io::Writer<Concurrency>>> file_ptrs_;
};

#ifdef DCA_HAVE_MPI
extern template class SignalHandler<dca::parallel::MPIConcurrency>;
#endif

}  // namespace util
}  // namespace dca

#endif  // DCA_UTIL_SIGNAL_HANDLER_HPP
