//-*-C++-*-

#ifndef DCA_QMCI_POSIX_JACKET_FOR_MC_ACCUMULATION_H
#define DCA_QMCI_POSIX_JACKET_FOR_MC_ACCUMULATION_H

#include <queue>
#include <iostream>
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/qmci_accumulator.h"

namespace DCA {
namespace QMCI {
/*!
 * \ingroup POSIX-TEMPLATES
 * \brief   A posix-jacket that implements a MC-accumulator, independent of the MC-TYPE.
 * \author  Peter Staar, Raffaele Solca
 * \version 1.0
 */
template <class qmci_accumulator_type>
class posix_qmci_accumulator : protected qmci_accumulator_type {
  typedef typename qmci_accumulator_type::my_parameters_type parameters_type;
  typedef typename qmci_accumulator_type::my_MOMS_type MOMS_type;

  typedef posix_qmci_accumulator<qmci_accumulator_type> this_type;

public:
  posix_qmci_accumulator(parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id);

  ~posix_qmci_accumulator();

  template <typename walker_type>
  void update_from(walker_type& walker);

  void wait_for_qmci_walker();

  void measure(pthread_mutex_t& mutex_queue, std::queue<this_type*>& accumulators_queue);

  // Sums all accumulated objects of this accumulator to the equivalent objects of the 'other'
  // accumulator.
  void sum_to(qmci_accumulator_type& other);

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

}  // namespace QMCI
}  // namespace DCA

#endif  // DCA_QMCI_POSIX_JACKET_FOR_MC_ACCUMULATION_H
