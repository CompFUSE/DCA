// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (peter.w.j.staar@gmail.com)
//         Raffaele Solca' (rasolca@itp.phys.ethz.ch)
//
// A posix jacket that implements a MC accumulator independent of the MC method.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_POSIX_QMCI_ACCUMULATOR_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_POSIX_QMCI_ACCUMULATOR_H

#include <pthread.h>
#include <queue>
#include <stdexcept>

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <class qmci_accumulator_type>
class posix_qmci_accumulator : protected qmci_accumulator_type {
  typedef typename qmci_accumulator_type::my_parameters_type parameters_type;
  typedef typename qmci_accumulator_type::my_MOMS_type MOMS_type;

  typedef posix_qmci_accumulator<qmci_accumulator_type> this_type;

public:
  posix_qmci_accumulator(parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id);

  ~posix_qmci_accumulator();

  using qmci_accumulator_type::initialize;
  using qmci_accumulator_type::finalize;
  // using qmci_accumulator_type::to_JSON;
  using qmci_accumulator_type::get_configuration;

  template <typename walker_type>
  void update_from(walker_type& walker);

  void wait_for_qmci_walker();

  void measure(pthread_mutex_t& mutex_queue, std::queue<this_type*>& accumulators_queue);

  // void sum_to(qmci_accumulator_type& accumulator_obj);
  // int get_expansion_order();

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sum_to(qmci_accumulator_type& other);

protected:
  using qmci_accumulator_type::get_Gflop;
  using qmci_accumulator_type::get_sign;
  using qmci_accumulator_type::get_number_of_measurements;

private:
  using qmci_accumulator_type::parameters;
  using qmci_accumulator_type::MOMS;

  /*
      using qmci_accumulator_type::HS_configuration_e_UP;
      using qmci_accumulator_type::HS_configuration_e_DN;

      using qmci_accumulator_type::visited_expansion_order_k;

      using qmci_accumulator_type::K_r_t;

      using qmci_accumulator_type::G_r_t;
      using qmci_accumulator_type::G_r_t_stddev;

      using qmci_accumulator_type::charge_cluster_moment;
      using qmci_accumulator_type::magnetic_cluster_moment;
      using qmci_accumulator_type::dwave_pp_correlator;

      using qmci_accumulator_type::M_r_w;
      using qmci_accumulator_type::M_r_w_squared;

      using qmci_accumulator_type::G4;
  */

  int thread_id;

  bool measuring;
  pthread_cond_t start_measuring;

  pthread_mutex_t mutex_accumulator;
};

template <class qmci_accumulator_type>
posix_qmci_accumulator<qmci_accumulator_type>::posix_qmci_accumulator(parameters_type& parameters_ref,
                                                                      MOMS_type& MOMS_ref, int id)
    : qmci_accumulator_type(parameters_ref, MOMS_ref, id), thread_id(id), measuring(false) {
  pthread_cond_init(&start_measuring, NULL);
  pthread_mutex_init(&mutex_accumulator, NULL);
}

template <class qmci_accumulator_type>
posix_qmci_accumulator<qmci_accumulator_type>::~posix_qmci_accumulator() {
  pthread_cond_destroy(&start_measuring);
  pthread_mutex_destroy(&mutex_accumulator);
}

//     template<class qmci_accumulator_type>
//     int posix_qmci_accumulator<qmci_accumulator_type>::get_expansion_order()
//     {
//       return HS_configuration_e_UP.size();
//     }

template <class qmci_accumulator_type>
template <typename walker_type>
void posix_qmci_accumulator<qmci_accumulator_type>::update_from(walker_type& walker) {
  {
    pthread_mutex_lock(&mutex_accumulator);

    if (measuring)
      throw std::logic_error(__FUNCTION__);

    qmci_accumulator_type::update_from(walker);

    measuring = true;

    pthread_cond_signal(&start_measuring);

    pthread_mutex_unlock(&mutex_accumulator);
  }
}

template <class qmci_accumulator_type>
void posix_qmci_accumulator<qmci_accumulator_type>::wait_for_qmci_walker() {
  pthread_mutex_lock(&mutex_accumulator);

  while (!measuring)
    pthread_cond_wait(&start_measuring, &mutex_accumulator);

  pthread_mutex_unlock(&mutex_accumulator);
}

template <class qmci_accumulator_type>
void posix_qmci_accumulator<qmci_accumulator_type>::measure(
    pthread_mutex_t& /*mutex_queue*/, std::queue<this_type*>& /*accumulators_queue*/) {
  pthread_mutex_lock(&mutex_accumulator);

  qmci_accumulator_type::measure();

  measuring = false;

  pthread_mutex_unlock(&mutex_accumulator);
}

template <class qmci_accumulator_type>
void posix_qmci_accumulator<qmci_accumulator_type>::sum_to(qmci_accumulator_type& other) {
  pthread_mutex_lock(&mutex_accumulator);
  qmci_accumulator_type::sum_to(other);
  pthread_mutex_unlock(&mutex_accumulator);
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_PTHREAD_JACKET_POSIX_QMCI_ACCUMULATOR_H
