// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides an utility to silence the lapack library warnings.

#ifndef DCA_LINALG_LAPACK_SILENCE_LAPACK_HPP
#define DCA_LINALG_LAPACK_SILENCE_LAPACK_HPP

namespace dca {
namespace linalg {
namespace lapack {
// dca::linalg::lapack::

void silenceLapack();

}  // lapack
}  // linalg
}  // dca

#endif  // DCA_LINALG_LAPACK_SILENCE_LAPACK_HPP
