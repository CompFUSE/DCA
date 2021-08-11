// Copyright (C) 2021 ETH Zurich
// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// This file defines the Monte Carlo walker device. It's acceptable to do
// this copy to "generic" header because we will never be running in a mixed gpu environment
// that's just not how HPC systems are built.

#ifndef DCA_CONFIG_WALKER_DEVICE_HPP
#define DCA_CONFIG_WALKER_DEVICE_HPP

#include "dca/linalg/device_type.hpp"
#include "dca/linalg/util/info_hip.hpp"    // Declares printInfoDevices().
#include "dca/linalg/util/util_hipblas.hpp"  // Declares initializeMagma().

constexpr dca::linalg::DeviceType walker_device = dca::linalg::GPU;

#endif  // DCA_CONFIG_WALKER_DEVICE_HPP
