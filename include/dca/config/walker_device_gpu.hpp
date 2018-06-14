// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file defines the Monte Carlo walker device.

#ifndef DCA_CONFIG_WALKER_DEVICE_HPP
#define DCA_CONFIG_WALKER_DEVICE_HPP

#include "dca/linalg/device_type.hpp"
#include "dca/linalg/util/info_cuda.hpp"    // Declares printInfoDevices().
#include "dca/linalg/util/util_cublas.hpp"  // Declares initializeMagma().

constexpr dca::linalg::DeviceType walker_device = dca::linalg::GPU;

#endif  // DCA_CONFIG_WALKER_DEVICE_HPP
