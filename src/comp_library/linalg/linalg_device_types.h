//-*-C++-*-

#ifndef LIN_ALG_DEVICE_TYPES
#define LIN_ALG_DEVICE_TYPES

#include "dca/linalg/device_type.hpp"

namespace LIN_ALG {

// TODO: remove device_type and replace with DeviceType
typedef dca::linalg::DeviceType device_type;
const device_type CPU = dca::linalg::CPU;
const device_type GPU = dca::linalg::GPU;

enum copy_concurrency_types { SYNCHRONOUS, ASYNCHRONOUS };
typedef copy_concurrency_types copy_concurrency_type;
}

#endif
