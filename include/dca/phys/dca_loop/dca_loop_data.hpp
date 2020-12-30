// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class contains physical and computational data of the DCA(+) loop.

#ifndef DCA_PHYS_DCA_LOOP_DCA_LOOP_DATA_HPP
#define DCA_PHYS_DCA_LOOP_DCA_LOOP_DATA_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/filesystem.hpp"
#include "dca/io/reader.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/dca_iteration_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"

namespace dca {
namespace phys {
// dca::phys::

template <typename ParametersType>
class DcaLoopData {
public:
  using Concurrency = typename ParametersType::concurrency_type;
  using expansion_dmn_t = func::dmn_0<func::dmn<32, int>>;
  using DCA_iteration_domain_type = func::dmn_0<domains::DCA_iteration_domain>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using k_DCA =
      func::dmn_0<domains::cluster_domain<double, ParametersType::lattice_type::DIMENSION, domains::CLUSTER,
                                          domains::MOMENTUM_SPACE, domains::BRILLOUIN_ZONE>>;

  DcaLoopData();

  template <typename Writer>
  void write(Writer& writer);

  // Attempts to read the loop functions from 'filename'. If successful returns
  // the last completed iteration from the input file, otherwise it returns -1.
  int tryToRead(const std::string& filename, const std::string& format,
                const Concurrency& concurrency);

  func::function<double, DCA_iteration_domain_type> Gflop_per_mpi_task;
  func::function<double, DCA_iteration_domain_type> times_per_mpi_task;

  func::function<double, DCA_iteration_domain_type> thermalization_per_mpi_task;
  func::function<double, DCA_iteration_domain_type> MC_integration_per_mpi_task;

  func::function<double, DCA_iteration_domain_type> Gflops_per_mpi_task;
  func::function<double, DCA_iteration_domain_type> max_Gflops_per_mpi_task;

  func::function<double, DCA_iteration_domain_type> sign;

  func::function<double, DCA_iteration_domain_type> L2_Sigma_difference;

  func::function<double, func::dmn_variadic<nu, k_DCA, DCA_iteration_domain_type>> Sigma_zero_moment;
  func::function<double, func::dmn_variadic<nu, k_DCA, DCA_iteration_domain_type>> standard_deviation;

  func::function<std::complex<double>, func::dmn_variadic<nu, nu, expansion_dmn_t, DCA_iteration_domain_type>>
      sigma_lambda;

  func::function<double, func::dmn_variadic<nu, k_DCA, DCA_iteration_domain_type>> n_k;
  func::function<double, func::dmn_variadic<nu, k_DCA, DCA_iteration_domain_type>> A_k;
  func::function<double, func::dmn_variadic<nu, DCA_iteration_domain_type>> orbital_occupancies;

  func::function<double, DCA_iteration_domain_type> density;
  func::function<double, DCA_iteration_domain_type> chemical_potential;
  func::function<double, DCA_iteration_domain_type> average_expansion_order;

  int last_completed_iteration = -1;
};

template <typename ParametersType>
DcaLoopData<ParametersType>::DcaLoopData()
    : Gflop_per_mpi_task("Gflop_per_mpi_task"),
      times_per_mpi_task("times_per_mpi_task"),

      thermalization_per_mpi_task("thermalization_per_mpi_task"),
      MC_integration_per_mpi_task("MC_integration_per_mpi_task"),

      Gflops_per_mpi_task("Gflops_per_mpi_task"),
      max_Gflops_per_mpi_task("max_Gflops_per_mpi_task"),

      sign("sign"),
      L2_Sigma_difference("L2_Sigma_difference"),

      Sigma_zero_moment("Sigma_zero_moment"),
      standard_deviation("standard_deviation"),

      sigma_lambda("sigma_lambda"),

      n_k("n_k"),
      A_k("A_k"),
      orbital_occupancies("orbital-occupancies"),

      density("density"),
      chemical_potential("chemical-potential"),
      average_expansion_order("expansion_order") {}

template <typename ParametersType>
template <typename Writer>
void DcaLoopData<ParametersType>::write(Writer& writer) {
  writer.open_group("DCA-loop-functions");

  {
    writer.execute(Gflop_per_mpi_task);
    writer.execute(times_per_mpi_task);
    writer.execute(Gflops_per_mpi_task);

    writer.execute(sign);

    writer.execute(L2_Sigma_difference);
    writer.execute(standard_deviation);

    writer.execute(chemical_potential);
    writer.execute(density);
    writer.execute(average_expansion_order);

    writer.execute(Sigma_zero_moment);

    writer.execute(n_k);
    writer.execute(A_k);
    writer.execute(orbital_occupancies);

    writer.execute("completed-iteration", last_completed_iteration);
  }

  writer.close_group();
}

template <typename ParametersType>
int DcaLoopData<ParametersType>::tryToRead(const std::string& filename, const std::string& format,
                                           const Concurrency& concurrency) {
  if (concurrency.id() == concurrency.first() && filesystem::exists(filename)) {
    io::Reader reader(format);

    reader.open_file(filename);
    reader.open_group("DCA-loop-functions");

    reader.execute("completed-iteration", last_completed_iteration);

    reader.execute(Gflop_per_mpi_task);
    reader.execute(times_per_mpi_task);
    reader.execute(Gflops_per_mpi_task);

    reader.execute(sign);

    reader.execute(L2_Sigma_difference);
    reader.execute(standard_deviation);

    reader.execute(chemical_potential);
    reader.execute(density);
    reader.execute(average_expansion_order);

    reader.execute(Sigma_zero_moment);

    reader.execute(n_k);
    reader.execute(A_k);
    reader.execute(orbital_occupancies);

    reader.open_group("DCA-loop-functions");
    reader.close_file();
  }

  concurrency.broadcast(last_completed_iteration);
  return last_completed_iteration;
}

}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_LOOP_DCA_LOOP_DATA_HPP
