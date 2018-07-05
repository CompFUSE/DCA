// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This file provides the type of devices used for linear algebra routines.

#ifndef DCA_LINALG_DEVICE_TYPE_HPP
#define DCA_LINALG_DEVICE_TYPE_HPP

namespace dca {
namespace linalg {
// dca::linalg::

enum DeviceType : int { CPU, GPU };

}  // linalg
}  // dca

#endif  // DCA_LINALG_DEVICE_TYPE_HPP
