// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
//
// Wrapper to an instance of ADIOS2Reader, HDF5Reader or JSONReader.

#include "dca/io/reader.hpp"

namespace dca::io {

template class Reader<dca::parallel::NoConcurrency>;
#ifdef DCA_HAVE_MPI
template class Reader<dca::parallel::MPIConcurrency>;
#endif

}  // namespace dca::io
