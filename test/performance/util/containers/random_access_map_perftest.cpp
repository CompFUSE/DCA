// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// RandomAccessMap performance test

#include "dca/util/containers/random_access_map.hpp"

#include <map>
#include <unordered_map>
#include <random>
#include <string>

#include <benchmark/benchmark.h>

#if __has_include(<absl/container/btree_map.h>)
#define HAVE_ABSL
#include <absl/container/btree_map.h>
#endif

const unsigned n_init = 10000;
const unsigned n_test = 10;
std::vector<int> keys;
std::vector<int> vals;

void init() {
  static bool initialized = false;
  if (initialized)
    return;
  initialized = true;

  for (int i = 0; i < n_init + n_test; ++i) {
    keys.push_back(i);
    vals.push_back(i);
  }

  std::random_shuffle(keys.begin(), keys.end());
  std::random_shuffle(vals.begin(), vals.end());
}

template <template <class, class> class Map>
void performSTLTest(benchmark::State& state) {
  init();
  Map<int, int> map;

  for (int i = 0; i < state.range(0); ++i)
    map.insert(std::make_pair(keys[i], vals[i]));

  for (auto _ : state) {
    for (int i = n_init; i < n_init + n_test; ++i)
      map.insert(std::make_pair(keys[i], vals[i]));
    for (int i = n_init; i < n_init + n_test; ++i)
      map.erase(keys[i]);
  }
}

static void BM_StdMap(benchmark::State& state) {
  performSTLTest<std::map>(state);
}
BENCHMARK(BM_StdMap)->Arg(100)->Arg(1000)->Arg(10000);

static void BM_StdUMap(benchmark::State& state) {
  performSTLTest<std::unordered_map>(state);
}
BENCHMARK(BM_StdUMap)->Arg(100)->Arg(1000)->Arg(10000);

#ifdef HAVE_ABSL
static void BM_AbslMap(benchmark::State& state) {
  performSTLTest<absl::btree_map>(state);
}
BENCHMARK(BM_AbslMap)->Arg(100)->Arg(1000)->Arg(10000);
#endif

static void BM_MyMap(benchmark::State& state) {
  init();
  dca::util::RandomAccessMap<int, int> map;

  for (int i = 0; i < state.range(0); ++i)
    map.insert(keys[i], vals[i]);

  for (auto _ : state) {
    for (int i = n_init; i < n_init + n_test; ++i)
      map.insert(keys[i], vals[i]);
    for (int i = n_init; i < n_init + n_test; ++i)
      map.erase(keys[i]);
  }
}
BENCHMARK(BM_MyMap)->Arg(100)->Arg(1000)->Arg(10000);
