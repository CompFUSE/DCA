// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// Interface to the kernels used by SpaceTransform2DGpu

#ifndef DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_KERNELS_INTERFACE
#define DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_KERNELS_INTERFACE
#ifdef DCA_HAVE_CUDA

#include <complex>
#include <cuda.h>

namespace dca {
namespace math {
namespace transform {
namespace details {
// dca::math::transform::details::

template <typename Real>
void rearrangeResult(const std::complex<Real>* in, int ldi, std::complex<Real>* out, int ldo, int nb,
                int nk, int nw, cudaStream_t stream);

}  // details
}  // transform
}  // math
}  // dca

#endif  // DCA_HAVE_CUDA
#endif  // DCA_MATH_FUNCTION_TRANSFORM_SPECIAL_TRANSFORMS_SPACE_KERNELS_INTERFACE
