//-*-C++-*-

#ifndef BASIC_CUDA_FUNCTIONS_CU_H
#define BASIC_CUDA_FUNCTIONS_CU_H

#include <iostream>
#include <cuda_runtime.h>

cudaDeviceProp& get_device_properties() {
  static bool init = false;
  static cudaDeviceProp prop;

  if (!init) {
    int count;

    cudaGetDeviceCount(&count);

    if (count < 1)
      throw std::logic_error(__FUNCTION__);

    cudaGetDeviceProperties(&prop, 0);

    init = true;
  }

  return prop;
}

int get_number_of_threads() {
  const static int N_th = get_device_properties().maxThreadsPerBlock;
  return N_th;
}

#endif
