// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the bubble in the particle-hole and particle-particle channel.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_COMPUTE_BARE_BUBBLE_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_COMPUTE_BARE_BUBBLE_H

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "dca/parallel/util/get_bounds.hpp"
#include "dca/parallel/util/threading_data.hpp"
#include "dca/util/print_time.hpp"
#include "comp_library/function_library/include_function_library.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"

namespace DCA {

namespace SERIES_EXPANSION {

enum channel { ph, pp };
using channel_type = channel;

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
class compute_bubble {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;
  using Threading = typename parameters_type::ThreadingType;

  using w = dmn_0<frequency_domain>;
  using w_VERTEX_BOSONIC = dmn_0<DCA::vertex_frequency_domain<DCA::EXTENDED_BOSONIC>>;

  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;
  using nu = dmn_variadic<b, s>;  // orbital-spin index
  using b_b = dmn_variadic<b, b>;

  using G_function_type =
      FUNC_LIB::function<std::complex<double>, dmn_variadic<nu, nu, k_dmn_t, w_dmn_t>>;
  using function_type =
      FUNC_LIB::function<std::complex<double>, dmn_variadic<b_b, b_b, k_dmn_t, w_VERTEX_BOSONIC>>;

public:
  compute_bubble(parameters_type& parameters_ref);

  function_type& get_function();

  void execute_on_lattice(G_function_type& G);
  void execute_on_cluster(G_function_type& G);

  void threaded_execute_on_cluster(G_function_type& G);

  template <typename Writer>
  void write(Writer& writer);

private:
  void execute_on_lattice_ph(G_function_type& S);
  void execute_on_cluster_ph(G_function_type& G);

  static void* threaded_execute_on_cluster_ph(void* data);

  void execute_on_lattice_pp(G_function_type& S);
  void execute_on_cluster_pp(G_function_type& G);

  static void* threaded_execute_on_cluster_pp(void* data);

private:
  parameters_type& parameters;
  concurrency_type& concurrency;

protected:
  FUNC_LIB::function<std::complex<double>, dmn_variadic<b_b, b_b, k_dmn_t, w_VERTEX_BOSONIC>> chi;

private:
  struct bubble_data {
    G_function_type* G_ptr;
    function_type* chi_ptr;

    concurrency_type* concurrency_ptr;
  };
};

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::compute_bubble(
    parameters_type& parameters_ref)
    : parameters(parameters_ref), concurrency(parameters_ref.get_concurrency()) {}

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
template <typename Writer>
void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::write(Writer& writer) {
  writer.execute(chi);
}

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
typename compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::function_type& compute_bubble<
    channel_value, parameters_type, k_dmn_t, w_dmn_t>::get_function() {
  return chi;
}

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::execute_on_cluster(
    G_function_type& G) {
  switch (channel_value) {
    case ph:
      execute_on_cluster_ph(G);
      break;

    case pp:
      execute_on_cluster_pp(G);
      break;

    default:
      throw std::logic_error(__FUNCTION__);
  }
}

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::threaded_execute_on_cluster(
    G_function_type& G) {
  if (concurrency.id() == 0)
    std::cout << "\n\n\t\t"
              << "threaded_execute_on_cluster compute-bubble"
              << "\n\n";

  profiler_type profiler("threaded_execute_on_cluster compute-bubble", "HTS", __LINE__);

  {
    int nr_threads = parameters.get_nr_HTS_threads();

    bubble_data args;

    args.G_ptr = &G;
    args.chi_ptr = &chi;

    args.concurrency_ptr = &concurrency;

    Threading pthreads;
    switch (channel_value) {
      case ph:
        pthreads.execute(nr_threads, threaded_execute_on_cluster_ph, (void*)&args);
        break;

      case pp:
        pthreads.execute(nr_threads, threaded_execute_on_cluster_pp, (void*)&args);
        break;

      default:
        throw std::logic_error(__FUNCTION__);
    }
  }

  concurrency.sum(chi);

  {
    double factor = -1. / (parameters.get_beta() * k_dmn_t::dmn_size());

    chi *= factor;
  }
}

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::execute_on_cluster_ph(
    G_function_type& G) {
  // cout << __FUNCTION__ << endl;
  if (concurrency.id() == 0)
    std::cout << "\n\n\t\t ph-buble \n\n" << std::endl;

  chi.get_name() = "ph-bubble";

  chi = 0.;

  assert(std::fabs(w_VERTEX_BOSONIC::get_elements()[w_VERTEX_BOSONIC::dmn_size() / 2]) < 1.e-6);

  for (int q_ind = 0; q_ind < k_dmn_t::dmn_size(); ++q_ind) {
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
      int k_plus_q = k_dmn_t::parameter_type::add(q_ind, k_ind);

      for (int nu_ind = 0; nu_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu_ind) {
        int nu_c = (nu_ind - w_VERTEX_BOSONIC::dmn_size() / 2);

        for (int w_ind = std::fabs(nu_c); w_ind < w_dmn_t::dmn_size() - std::fabs(nu_c); ++w_ind) {
          int w_plus_nu = w_ind + nu_c;

          for (int j1 = 0; j1 < b::dmn_size(); ++j1)
            for (int j0 = 0; j0 < b::dmn_size(); ++j0)
              for (int i1 = 0; i1 < b::dmn_size(); ++i1)
                for (int i0 = 0; i0 < b::dmn_size(); ++i0)
                  chi(i0, i1, j0, j1, q_ind, nu_ind) +=
                      G(i0, j1, k_ind, w_ind) * G(i1, j0, k_plus_q, w_plus_nu);
        }
      }
    }
  }

  double factor = -1. / (parameters.get_beta() * k_dmn_t::dmn_size());
  chi *= factor;
}

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
void* compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::threaded_execute_on_cluster_ph(
    void* void_ptr) {
  dca::parallel::ThreadingData* data_ptr = static_cast<dca::parallel::ThreadingData*>(void_ptr);
  bubble_data* bubble_ptr = static_cast<bubble_data*>(data_ptr->arg);

  G_function_type& G = *(bubble_ptr->G_ptr);
  function_type& chi = *(bubble_ptr->chi_ptr);

  concurrency_type& concurrency = *(bubble_ptr->concurrency_ptr);

  k_dmn_t k_dmn;
  std::pair<int, int> k_bounds = concurrency.get_bounds(k_dmn);

  int id = data_ptr->id;
  int nr_threads = data_ptr->num_threads;

  k_dmn_t q_dmn;
  std::pair<int, int> q_bounds = dca::parallel::util::getBounds(id, nr_threads, q_dmn);

  for (int q_ind = q_bounds.first; q_ind < q_bounds.second; ++q_ind) {
    double percentage = double(q_ind - q_bounds.first) / double(q_bounds.second - q_bounds.first);

    if (concurrency.id() == 0 and id == 0 and (int(100 * percentage) % 10 == 0))
      std::cout << "\t" << int(100 * percentage) << " % finished\t" << dca::util::print_time()
                << "\n";

    for (int k_ind = k_bounds.first; k_ind < k_bounds.second; ++k_ind) {
      int k_plus_q = k_dmn_t::parameter_type::add(q_ind, k_ind);

      for (int nu_ind = 0; nu_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu_ind) {
        int nu_c = (nu_ind - w_VERTEX_BOSONIC::dmn_size() / 2);

        for (int w_ind = std::fabs(nu_c); w_ind < w_dmn_t::dmn_size() - std::fabs(nu_c); ++w_ind) {
          int w_plus_nu = w_ind + nu_c;

          for (int j1 = 0; j1 < b::dmn_size(); ++j1)
            for (int j0 = 0; j0 < b::dmn_size(); ++j0)
              for (int i1 = 0; i1 < b::dmn_size(); ++i1)
                for (int i0 = 0; i0 < b::dmn_size(); ++i0)
                  chi(i0, i1, j0, j1, q_ind, nu_ind) +=
                      G(i0, j1, k_ind, w_ind) * G(i1, j0, k_plus_q, w_plus_nu);
        }
      }
    }
  }

  return 0;
}

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
void compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::execute_on_cluster_pp(
    G_function_type& G) {
  // cout << __FUNCTION__ << endl;
  if (concurrency.id() == 0)
    std::cout << "\n\n\t\t pp-buble \n\n" << std::endl;

  chi.get_name() = "pp-bubble";

  chi = 0.;

  assert(std::fabs(w_VERTEX_BOSONIC::get_elements()[w_VERTEX_BOSONIC::dmn_size() / 2]) < 1.e-6);

  for (int q_ind = 0; q_ind < k_dmn_t::dmn_size(); ++q_ind) {
    for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); ++k_ind) {
      int q_minus_k = k_dmn_t::parameter_type::subtract(k_ind, q_ind);

      for (int nu_ind = 0; nu_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu_ind) {
        int nu_c = (nu_ind - w_VERTEX_BOSONIC::dmn_size() / 2);

        for (int w_ind = std::fabs(nu_c); w_ind < w_dmn_t::dmn_size() - std::fabs(nu_c); ++w_ind) {
          int nu_minus_w = nu_c + (w::dmn_size() - 1 - w_ind);

          for (int j1 = 0; j1 < b::dmn_size(); ++j1)
            for (int j0 = 0; j0 < b::dmn_size(); ++j0)
              for (int i1 = 0; i1 < b::dmn_size(); ++i1)
                for (int i0 = 0; i0 < b::dmn_size(); ++i0)
                  chi(i0, i1, j0, j1, q_ind, nu_ind) +=
                      G(i0, j0, k_ind, w_ind) * G(i1, j1, q_minus_k, nu_minus_w);
        }
      }
    }
  }

  double factor = -1. / (parameters.get_beta() * k_dmn_t::dmn_size());
  chi *= factor;
}

template <channel_type channel_value, class parameters_type, class k_dmn_t, class w_dmn_t>
void* compute_bubble<channel_value, parameters_type, k_dmn_t, w_dmn_t>::threaded_execute_on_cluster_pp(
    void* void_ptr) {
  dca::parallel::ThreadingData* data_ptr = static_cast<dca::parallel::ThreadingData*>(void_ptr);
  bubble_data* bubble_ptr = static_cast<bubble_data*>(data_ptr->arg);

  G_function_type& G = *(bubble_ptr->G_ptr);
  function_type& chi = *(bubble_ptr->chi_ptr);

  concurrency_type& concurrency = *(bubble_ptr->concurrency_ptr);

  k_dmn_t k_dmn;
  std::pair<int, int> k_bounds = concurrency.get_bounds(k_dmn);

  int id = data_ptr->id;
  int nr_threads = data_ptr->num_threads;

  k_dmn_t q_dmn;
  std::pair<int, int> q_bounds = dca::parallel::util::getBounds(id, nr_threads, q_dmn);

  for (int q_ind = q_bounds.first; q_ind < q_bounds.second; ++q_ind) {
    double percentage = double(q_ind - q_bounds.first) / double(q_bounds.second - q_bounds.first);

    if (concurrency.id() == 0 and id == 0 and (int(100 * percentage) % 10 == 0))
      std::cout << "\t" << int(100 * percentage) << " % finished\t" << dca::util::print_time()
                << "\n";

    for (int k_ind = k_bounds.first; k_ind < k_bounds.second; ++k_ind) {
      int q_minus_k = k_dmn_t::parameter_type::subtract(k_ind, q_ind);

      for (int nu_ind = 0; nu_ind < w_VERTEX_BOSONIC::dmn_size(); ++nu_ind) {
        int nu_c = (nu_ind - w_VERTEX_BOSONIC::dmn_size() / 2);

        for (int w_ind = std::fabs(nu_c); w_ind < w_dmn_t::dmn_size() - std::fabs(nu_c); ++w_ind) {
          int nu_minus_w = nu_c + (w::dmn_size() - 1 - w_ind);

          for (int j1 = 0; j1 < b::dmn_size(); ++j1)
            for (int j0 = 0; j0 < b::dmn_size(); ++j0)
              for (int i1 = 0; i1 < b::dmn_size(); ++i1)
                for (int i0 = 0; i0 < b::dmn_size(); ++i0)
                  chi(i0, i1, j0, j1, q_ind, nu_ind) +=
                      G(i0, j0, k_ind, w_ind) * G(i1, j1, q_minus_k, nu_minus_w);
        }
      }
    }
  }

  return 0;
}
}
}

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_SERIES_EXPANSION_COMPUTE_BARE_BUBBLE_H
