// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides a wrapper class for the random number library of C++11.
// The class is templated on the random number engine and can be instantiated for example with
// std::mt19937_64 or std::ranlux48.
// It produces uniformly distributed pseudo-random numbers of type double in the interval [0, 1).
// In addition, this class keeps track of the number of created objects in order to assign a unique
// ID to each of them, that is then used to generate a unique seed.

#ifndef DCA_MATH_RANDOM_STD_RANDOM_WRAPPER_HPP
#define DCA_MATH_RANDOM_STD_RANDOM_WRAPPER_HPP

#include <cstdint>  // for uint64_t
#include <random>
#include "dca/math/random/random_utils.hpp"

namespace dca {
namespace math {
namespace random {
// dca::math::random::

template <typename Engine>
class StdRandomWrapper {
public:
  StdRandomWrapper(const int proc_id, const int num_procs, const uint64_t seed = 0)
      : local_id_(counter_++),
        global_id_(detail::getGlobalId(local_id_, proc_id, num_procs)),
        initial_seed_(seed),
        seed_(detail::generateSeed(global_id_, seed)),
        engine_(seed_),
        distro_(0., 1.) {}

  // Make the random number generator object non-copyable, but move-constructible.
  // The implicit move assignment operator is deleted since with we have non-static const members.
  StdRandomWrapper(const StdRandomWrapper&) = delete;
  StdRandomWrapper& operator=(const StdRandomWrapper&) = delete;
  StdRandomWrapper(StdRandomWrapper&&) = default;
  // StdRandomWrapper& operator=(StdRandomWrapper&&) = delete;

  ~StdRandomWrapper() = default;

  inline int getGlobalId() const {
    return global_id_;
  }

  inline uint64_t getInitialSeed() const {
    return initial_seed_;
  }

  inline uint64_t getSeed() const {
    return seed_;
  }

  // Reset the static counter. For testing purposes.
  static void resetCounter() {
    counter_ = 0;
  }

  // Returns a uniformly distributied pseudo-random number in the interval [0, 1).
  inline double operator()() {
    return distro_(engine_);
  }

private:
  static int counter_;

  const int local_id_;
  const int global_id_;
  const uint64_t initial_seed_;
  const uint64_t seed_;

  Engine engine_;
  std::uniform_real_distribution<double> distro_;
};

template <typename Engine>
int StdRandomWrapper<Engine>::counter_ = 0;

// Partial specialization for std::linear_congruential_engine in order to prevent usage.
template <typename UIntType, UIntType a, UIntType c, UIntType m>
class StdRandomWrapper<std::linear_congruential_engine<UIntType, a, c, m>>;

}  // random
}  // math
}  // dca

#endif  // DCA_MATH_RANDOM_STD_RANDOM_WRAPPER_HPP
