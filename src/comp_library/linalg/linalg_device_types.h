//-*-C++-*-

#ifndef LIN_ALG_DEVICE_TYPES
#define LIN_ALG_DEVICE_TYPES

namespace LIN_ALG {

  enum    device_types {CPU, GPU};
  typedef device_types device_type;

  enum    copy_concurrency_types {SYNCHRONOUS, ASYNCHRONOUS};
  typedef copy_concurrency_types copy_concurrency_type;
}

#endif 
