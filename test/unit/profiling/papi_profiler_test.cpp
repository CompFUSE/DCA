#include "dca/profiling/counting_profiler.hpp"
#include "dca/profiling/events/papi_and_time_event.hpp"

#include <vector>
#include <future>

#include "gtest/gtest.h"

#include "dca/io/json/json_reader.hpp"

using Profiler = dca::profiling::CountingProfiler<dca::profiling::papi_and_time_event<long long int>>;

TEST(MyPapiProfilerTest, Parallel) {
  Profiler::start();
  constexpr int n_threads = 4;
  constexpr int n = 1e3;

  {
    std::vector<std::future<void>> futures;
    for (int id = 0; id < n_threads; ++id)
      futures.emplace_back(std::async(std::launch::async, [id]() {
        Profiler::start_threading(id);

        std::vector<float> a(n, 1), b(n, 2), c(n, 1);

        {
          Profiler prof(__FUNCTION__, "MyPapiProfilerTest", __LINE__, id);
          for (int i = 0; i < n; ++i)
            c[i] = a[i] + b[i];
        }

        Profiler::stop_threading(id);
      }));
  }

  Profiler::stop("profile.json");

  // Read the output and check the number of flops.
  dca::io::JSONReader reader;
  reader.open_file("profile.json");
  reader.open_group("0");

  int flops;
  reader.execute("PAPI_FP_OPS", flops);

  EXPECT_EQ(n_threads * n, flops);
}
