//-*-C++-*-                                                                                                                                                                                                                                                                         

#ifndef LIN_ALG_CUBLAS_MANAGER_TEM_H
#define LIN_ALG_CUBLAS_MANAGER_TEM_H

namespace LIN_ALG {

  template<device_type device_name_0, device_type device_name_1>
  struct CUBLAS_DEVICE_NAME
  {};

  template<>
  struct CUBLAS_DEVICE_NAME<CPU, CPU>
  {
    const static device_type device_t = CPU;
  };

  template<>
  struct CUBLAS_DEVICE_NAME<GPU, CPU>
  {
    const static device_type device_t = GPU;
  };

  template<>
  struct CUBLAS_DEVICE_NAME<CPU, GPU>
  {
    const static device_type device_t = GPU;
  };

  template<>
  struct CUBLAS_DEVICE_NAME<GPU, GPU>
  {
    const static device_type device_t = GPU;
  };
 

  template<device_type device_name>
  class CUBLAS_THREAD_MANAGER
  {};
  
}

#endif
