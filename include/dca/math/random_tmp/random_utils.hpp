// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides utility functions for the creation of seeds and random numbers.

#ifndef DCA_MATH_RANDOM_RANDOM_UTILS_HPP
#define DCA_MATH_RANDOM_RANDOM_UTILS_HPP

#include <cstdint>  // for uint64_t

namespace dca {
namespace math {
namespace random {
namespace detail {
// dca::math::random::detail::

// Computes from the process ID and the local ID (=ID within process) a unique global ID.
int getGlobalId(const int local_id, const int process_id, const int num_procs);

// Computes a unique seed from the global ID.
// The parameter 'offset' can be used to generate different seeds.
uint64_t generateSeed(const int global_rng_id, const uint64_t offset = 0);

// Hash function for 64 bit integer.
// Based on MurmurHash3's 64-bit Finalizer
// (https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp),
// with optimized parameters 'Mix13' from
// http://zimbry.blogspot.ch/2011/09/better-bit-mixing-improving-on.html.
uint64_t hash(uint64_t key);

}  // detail
}  // random
}  // math
}  // dca

#endif  // DCA_MATH_RANDOM_RANDOM_UTILS_HPP
