// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
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
