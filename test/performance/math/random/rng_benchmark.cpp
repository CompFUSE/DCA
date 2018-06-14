// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file compares the performance (in terms of speed) of different random number generators.

#include <iostream>
#include <numeric>  // for std::accumulate

#include "dca/math/random/random.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/util/print_type.hpp"

template <typename Generator>
double runBenchmark(const int draws, const int reps) {
  std::vector<double> timings;

  // Sum up all random numbers and print the result, otherwise the compiler optimizes everything
  // away.
  double sum = 0;

  for (int r = 0; r < reps; ++r) {
    dca::profiling::WallTime start;

    Generator rng(0, 1);
    for (int i = 0; i < draws; ++i)
      sum += rng();

    dca::profiling::WallTime end;
    dca::profiling::Duration time(end, start);

    timings.push_back(time.sec + 1.e-6 * time.usec);
  }

  double avg_time = std::accumulate(timings.begin(), timings.end(), 0.) / timings.size();

  std::cout << dca::util::Type<Generator>::print() << "\t" << avg_time << std::endl;

  return sum;
}

int main() {
  const unsigned int draws = 1e8;
  const int reps = 3;

  std::cout << "Running benchmark for " << draws << " draws with " << reps << " repetitions ...\n"
            << std::endl;
  std::cout << "Generator\t\t\taverage time [s]" << std::endl;

  // Standard random number library
  runBenchmark<dca::math::random::StdRandomWrapper<std::mt19937>>(draws, reps);
  runBenchmark<dca::math::random::StdRandomWrapper<std::mt19937_64>>(draws, reps);
  runBenchmark<dca::math::random::StdRandomWrapper<std::ranlux48_base>>(draws, reps);
  runBenchmark<dca::math::random::StdRandomWrapper<std::ranlux48>>(draws, reps);  // very slow!

  return 0;
}
