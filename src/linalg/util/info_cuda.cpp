// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// This file implements info_cuda functions.

#include "dca/linalg/util/error_cuda.hpp"
#include <cuda_runtime.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>

namespace dca {
namespace linalg {
namespace util {
// dca::linalg::util::

void printInfoDevices() {
  int nr_devices;
  cudaGetDeviceCount(&nr_devices);

  std::stringstream s;
  s << "\n"
    << "********************************************************************************\n"
    << "**********                            CUDA                            **********\n"
    << "********************************************************************************\n"
    << "\n"
    << "CUDA found " << nr_devices << " devices."
    << "\n";

  for (int i = 0; i < nr_devices; ++i) {
    s << "\n";
    s << "  Device " << i << ":"
      << "\n";

    cudaDeviceProp dev_prop;
    cudaGetDeviceProperties(&dev_prop, i);

    s << "  Name:                      " << dev_prop.name << "\n";
    s << "  Compute capability:        " << dev_prop.major << "." << dev_prop.minor << "\n";
    s << "  Number of multiprocessors: " << dev_prop.multiProcessorCount << "\n";
    s << "  Global memory:             " << dev_prop.totalGlobalMem << " bytes"
      << "\n";
    s << "  Constant memory:           " << dev_prop.totalConstMem << " bytes"
      << "\n";
    s << "  Shared memory per block:   " << dev_prop.sharedMemPerBlock << " bytes"
      << "\n";
    s << "  Registers per block:       " << dev_prop.regsPerBlock << "\n";
    s << "  Maximum memory pitch:      " << dev_prop.memPitch << " bytes"
      << "\n";
    s << "  Warp size:                 " << dev_prop.warpSize << "\n";
    s << "  Maximum threads per block: " << dev_prop.maxThreadsPerBlock << "\n";

    s << "  Maximum size of blocks:    " << dev_prop.maxThreadsDim[0] << " x "
      << dev_prop.maxThreadsDim[1] << " x " << dev_prop.maxThreadsDim[2] << "\n";

    s << "  Maximum size of grids:     " << dev_prop.maxGridSize[0] << " x "
      << dev_prop.maxGridSize[1] << " x " << dev_prop.maxGridSize[2] << "\n";

    s << "  Clock frequency:           " << dev_prop.clockRate << " KHz"
      << "\n";
    s << "  Async engine count:        " << dev_prop.asyncEngineCount << "\n";
    s << "  Kernel execution timeout:  " << (dev_prop.kernelExecTimeoutEnabled ? "Yes" : "No")
      << "\n";
  }

  std::cout << s.str() << std::endl;
}

}  // util
}  // linalg
}  // dca
