// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file provides specific tests for the ADIOS2 reader and writer.

#include "dca/io/adios2/adios2_reader.hpp"
#include "dca/io/adios2/adios2_writer.hpp"
#include "dca/parallel/mpi_concurrency/mpi_concurrency.hpp"

#include <array>
#include <complex>
#include <string>
#include <vector>
#include <typeinfo>  // operator typeid
#include <iostream>
#include <sstream>
#include <fstream>

#include "dca/parallel/util/get_workload.hpp"
#include "dca/profiling/events/time.hpp"
#include "dca/profiling/events/time_event.hpp"

template <class T>
std::string VectorToString(const std::vector<T>& v) {
  std::stringstream ss;
  ss << "[";
  for (size_t i = 0; i < v.size(); ++i) {
    if (i != 0)
      ss << ",";
    ss << v[i];
  }
  ss << "]";
  return ss.str();
}

struct CommEnvironment {
  int rank;
  int comm_size;
  dca::parallel::MPIConcurrency* concurrency_ptr;
  adios2::ADIOS& adios;
};

template <int RANKS, int DMN3MULT, typename ST>
int functionReadWrite(CommEnvironment comm_env) {
  constexpr int dmn3_mult = DMN3MULT;
  constexpr int dmn3 = dmn3_mult * RANKS;
  constexpr int dmn1 = 10;
  constexpr int dmn2 = 16;
  using Dmn1 = dca::func::dmn_0<dca::func::dmn<dmn1>>;
  using Dmn2 = dca::func::dmn_0<dca::func::dmn<dmn2>>;
  using Dmn3 = dca::func::dmn_0<dca::func::dmn<dmn3>>;
  using Dmn = dca::func::dmn_variadic<Dmn1, Dmn2, Dmn3>;
  using Scalar = ST;
  const std::string typeStr = typeid(ST).name();

  constexpr uint64_t element_count = dmn1 * dmn2 * dmn3;
  constexpr uint64_t kbytes = element_count * sizeof(ST) / 1024;
  
  if (comm_env.comm_size != RANKS) {
    std::cerr << __func__ << " fast index not divisible by num ranks, "
              << std::to_string(comm_env.comm_size)
              << ", is not a good number of ranks to test this.\n"
              << "ADIOS2 would have UB." << std::endl;
    return -1;
  }

  dca::func::function<Scalar, Dmn> f1("parallelFunc");
  size_t dmn_size = 1;
  for (int l = 0; l < f1.signature(); ++l)
    dmn_size *= f1[l];

  int val = comm_env.rank * dmn_size;

  uint64_t start = 0;
  uint64_t end = 0;
  // This returns the linearized bounds of the function for a rank.
  dca::parallel::util::getComputeRange(comm_env.concurrency_ptr->id(),
                                       comm_env.concurrency_ptr->number_of_processors(), static_cast<uint64_t>(f1.size()),
                                       start, end);

  // only set this ranks values
  for (int i = start; i <= end; ++i)
    f1.data()[i] = ++val;

  // get the N-dimensional decomposition
  std::vector<int> subind_start = f1.linind_2_subind(start);
  std::vector<int> subind_end = f1.linind_2_subind(end);
  std::cout << "rank: " << comm_env.rank << " start lin = " << start
            << " subindicies = " << VectorToString(subind_start) << " end lin = " << end
            << " subindicies = " << VectorToString(subind_end) << std::endl;

  dca::profiling::WallTime start_write;

  const std::string fname("ADIOS2ParallelIOTest_" + typeStr + ".bp");
  {
    dca::io::ADIOS2Writer writer(comm_env.adios, comm_env.concurrency_ptr, true);
    writer.open_file(fname, true);

    // Because the caller needs to know if its function is distributed or not we will assume this is
    // so for the API as well. in the future I think something more sophisticated needs to be done
    // and the function will need to know its distribution, but for now we distribute only over the
    // fastest index

    writer.execute(f1, subind_start, subind_end);

    writer.close_file();
  }
  dca::profiling::WallTime end_write;

  dca::profiling::WallTime start_read;
  {
    // Read test file.
    if (!comm_env.rank) {
      std::cout << " Read back data with 3D selection " << std::endl;
    }

    dca::io::ADIOS2Reader reader(comm_env.adios, comm_env.concurrency_ptr, true);
    reader.open_file(fname);

    dca::func::function<Scalar, Dmn> f2("parallelFunc");

    if (!comm_env.rank) {
      std::cout << " Read back data with linear 1D selection " << std::endl;
    }

    /* TODO: This should be working on every rank */

    reader.close_file();
  }
  dca::profiling::WallTime end_read;

  std::ostringstream perf_file_name;
  perf_file_name << "adios2_perf_" << dmn1 << "x" << dmn2 << "x" << dmn3 << "_over_"
                 << comm_env.comm_size << ".txt";

  std::ofstream perf_file(perf_file_name.str());

  auto duration = [](dca::profiling::WallTime end, dca::profiling::WallTime start) {
    dca::profiling::Duration elapsed(end, start);
    return elapsed.sec + 1e-6 * elapsed.usec;
  };

  auto report = [&](std::ostream& o) -> std::ostream& {
    o << "Adios2 Performance on dca::function " << dmn1 << "x" << dmn2 << "x" << dmn3 << " over "
      << comm_env.comm_size << "ranks\n";
    auto write_duration = duration(end_write, start_write);
    auto write_speed = kbytes / write_duration;
    o << "write: " << write_duration << " s  " << write_speed / 1024 << " mb/s\n";
    auto read_duration = duration(end_read, start_read);
    auto read_speed = kbytes / read_duration;
    o << "read: " << read_duration << " s  " << read_speed / 1024 << " mb/s\n";
    return o;
  };
  report(std::cout);
  report(perf_file);

  return 0;
}

int main(int argc, char** argv) {
  dca::parallel::MPIConcurrency concurrency(argc, argv);
  adios2::ADIOS adios("", concurrency.get());

  CommEnvironment comm_env{concurrency.id(), concurrency.number_of_processors(), &concurrency, adios};
  
  functionReadWrite<ADIOS2_PERF_TEST_RANKS, 10, double>(comm_env);
  functionReadWrite<ADIOS2_PERF_TEST_RANKS, 100, double>(comm_env);
  functionReadWrite<ADIOS2_PERF_TEST_RANKS, 1000, double>(comm_env);
  functionReadWrite<ADIOS2_PERF_TEST_RANKS, 10000, double>(comm_env);
  functionReadWrite<ADIOS2_PERF_TEST_RANKS, 100000, double>(comm_env);

}
