//-*-C++-*-

#ifndef DCA_QMCI_POSIX_MC_INTEGRATOR_FOR_MC_H
#define DCA_QMCI_POSIX_MC_INTEGRATOR_FOR_MC_H

#include "dca/math_library/random_number_library//random_number_library.hpp"
#include "rng_type.inc"

namespace DCA {
/*!
 *  \defgroup POSIX-TEMPLATES
 */

/*!
 * \ingroup POSIX-TEMPLATES
 * \brief   A posix-MC_integrator that implements a threaded MC-integration, independent of the
 * MC-TYPE.
 * \author  Raffaele Solca, Peter Staar
 * \version 1.0
 */
template <class qmci_integrator_type>
class posix_qmci_integrator : protected qmci_integrator_type {
  typedef typename qmci_integrator_type::this_MOMS_type MOMS_type;
  typedef typename qmci_integrator_type::this_parameters_type parameters_type;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  using rng_type = random_number_generator;

  typedef typename qmci_integrator_type::walker_type walker_type;
  typedef typename qmci_integrator_type::accumulator_type accumulator_type;

  typedef posix_qmci_integrator<qmci_integrator_type> this_type;
  typedef QMCI::posix_qmci_accumulator<accumulator_type> posix_accumulator_type;

public:
  posix_qmci_integrator(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  ~posix_qmci_integrator();

  template <IO::FORMAT DATA_FORMAT>
  void write(IO::writer<DATA_FORMAT>& reader);

  void initialize(int dca_iteration0);

  void integrate();

  template <typename dca_info_struct_t>
  double finalize(dca_info_struct_t& dca_info_struct);

private:
  void print_thread_distribution();

  void set_the_rngs();

  static void* start_walker_static(void* arg);
  static void* start_accumulator_static(void* arg);

  void start_walker(int id);
  void start_accumulator(int id);

  void warm_up(walker_type& walker, int id);

  using qmci_integrator_type::compute_error_bars;
  using qmci_integrator_type::sum_measurements;
  using qmci_integrator_type::symmetrize_measurements;

private:
  using qmci_integrator_type::parameters;
  using qmci_integrator_type::MOMS;
  using qmci_integrator_type::concurrency;

  using qmci_integrator_type::total_time;

  using qmci_integrator_type::DCA_iteration;

  using qmci_integrator_type::accumulator;

  int dca_iteration;

  int acc_finished;

  int nr_walkers;
  int nr_accumulators;

  std::vector<rng_type> rng_vector;

  std::queue<posix_accumulator_type*> accumulators_queue;

  pthread_mutex_t mutex_print;
  pthread_mutex_t mutex_merge;
  pthread_mutex_t mutex_queue;
  pthread_mutex_t mutex_acc_finished;

  pthread_mutex_t mutex_numerical_error;
};

template <class qmci_integrator_type>
posix_qmci_integrator<qmci_integrator_type>::posix_qmci_integrator(parameters_type& parameters_ref,
                                                                   MOMS_type& MOMS_ref)
    : qmci_integrator_type(parameters_ref, MOMS_ref, false /*do not set rng*/),

      nr_walkers(parameters.get_nr_walkers()),
      nr_accumulators(parameters.get_nr_accumulators()),

      rng_vector(nr_walkers),
      accumulators_queue() {
  if (nr_walkers < 1 || nr_accumulators < 1) {
    std::cout << "\n\n\n"
              << "\t nr_walkers      --> " << nr_walkers << "\t nr_accumulators --> "
              << nr_accumulators << "\n\n\n"
              << "\t\t\t\t WTF are you doing !!!!\n\n";
    throw std::logic_error(__PRETTY_FUNCTION__);
  }

  // initialize the seeds and the rng's
  set_the_rngs();
}

template <class qmci_integrator_type>
posix_qmci_integrator<qmci_integrator_type>::~posix_qmci_integrator() {}

template <class qmci_integrator_type>
template <IO::FORMAT DATA_FORMAT>
void posix_qmci_integrator<qmci_integrator_type>::write(IO::writer<DATA_FORMAT>& writer) {
  qmci_integrator_type::write(writer);
}

template <class qmci_integrator_type>
void posix_qmci_integrator<qmci_integrator_type>::set_the_rngs() {
  int id_base = nr_walkers * concurrency.id();
  int tot_walkers = nr_walkers * concurrency.number_of_processors();
  for (int i = 0; i < nr_walkers; ++i)
    rng_vector[i].init_from_id(id_base + i, tot_walkers);
}

template <class qmci_integrator_type>
void posix_qmci_integrator<qmci_integrator_type>::initialize(int dca_iteration0) {
  dca_iteration = dca_iteration0;

  profiler_type profiler(__FUNCTION__, "posix-MC-Integration", __LINE__);

  qmci_integrator_type::initialize(dca_iteration);

  acc_finished = 0;

  pthread_mutex_init(&mutex_print, NULL);
  pthread_mutex_init(&mutex_merge, NULL);
  pthread_mutex_init(&mutex_queue, NULL);

  pthread_mutex_init(&mutex_acc_finished, NULL);
  pthread_mutex_init(&mutex_numerical_error, NULL);
}

template <class qmci_integrator_type>
void posix_qmci_integrator<qmci_integrator_type>::integrate() {
  profiler_type profiler(__FUNCTION__, "posix-MC-Integration", __LINE__);

  concurrency << "\n\t\t threaded QMC integration starts\n\n";

  std::vector<pthread_t> threads(nr_accumulators + nr_walkers);
  std::vector<std::pair<this_type*, int>> data(nr_accumulators + nr_walkers);

  {
    print_thread_distribution();

    PROFILER::WallTime start_time;

    int min_nr_walkers_nr_accumulators = std::min(nr_walkers, nr_accumulators);
    for (int i = 0; i < 2 * min_nr_walkers_nr_accumulators; ++i) {
      data[i] = std::pair<this_type*, int>(this, i);

      if (i % 2 == 0)
        pthread_create(&threads[i], NULL, start_walker_static, &data[i]);
      else
        pthread_create(&threads[i], NULL, start_accumulator_static, &data[i]);
    }

    for (int i = 2 * min_nr_walkers_nr_accumulators; i < nr_walkers + nr_accumulators; ++i) {
      data[i] = std::pair<this_type*, int>(this, i);

      if (min_nr_walkers_nr_accumulators != nr_walkers)
        pthread_create(&threads[i], NULL, start_walker_static, &data[i]);
      else
        pthread_create(&threads[i], NULL, start_accumulator_static, &data[i]);
    }

    void* rc;
    for (int i = 0; i < nr_walkers + nr_accumulators; ++i) {
      pthread_join(threads[i], &rc);
    }

    PROFILER::WallTime end_time;

    PROFILER::Duration duration(end_time, start_time);
    total_time = duration.sec + 1.e-6 * duration.usec;
  }

  symmetrize_measurements();

  compute_error_bars(parameters.get_number_of_measurements() * nr_accumulators);

  sum_measurements(parameters.get_number_of_measurements() * nr_accumulators);

  concurrency << "\n\t\t threaded QMC integration ends\n\n";
}

template <class qmci_integrator_type>
void posix_qmci_integrator<qmci_integrator_type>::print_thread_distribution() {
  if (concurrency.id() == 0) {
    int min_nr_walkers_nr_accumulators = std::min(nr_walkers, nr_accumulators);
    for (int i = 0; i < 2 * min_nr_walkers_nr_accumulators; ++i) {
      if (i % 2 == 0)
        std::cout << "\t pthread-id : " << i << "  -->   (walker)\n";
      else
        std::cout << "\t pthread-id : " << i << "  -->   (accumulator)\n";
    }

    for (int i = 2 * min_nr_walkers_nr_accumulators; i < nr_walkers + nr_accumulators; ++i) {
      if (min_nr_walkers_nr_accumulators != nr_walkers)
        std::cout << "\t pthread-id : " << i << "  -->   (walker)\n";
      else
        std::cout << "\t pthread-id : " << i << "  -->   (accumulator)\n";
    }
  }
}

template <class qmci_integrator_type>
template <typename dca_info_struct_t>
double posix_qmci_integrator<qmci_integrator_type>::finalize(dca_info_struct_t& dca_info_struct) {
  profiler_type profiler(__FUNCTION__, "posix-MC-Integration", __LINE__);

  double L2_Sigma_difference = qmci_integrator_type::finalize(dca_info_struct);

  pthread_mutex_destroy(&mutex_print);
  pthread_mutex_destroy(&mutex_merge);
  pthread_mutex_destroy(&mutex_queue);

  pthread_mutex_destroy(&mutex_acc_finished);
  pthread_mutex_destroy(&mutex_numerical_error);

  return L2_Sigma_difference;
}

template <class qmci_integrator_type>
void* posix_qmci_integrator<qmci_integrator_type>::start_walker_static(void* arg) {
  std::pair<this_type*, int>* data = reinterpret_cast<std::pair<this_type*, int>*>(arg);

  profiler_type::start_pthreading(data->second);

  data->first->start_walker(data->second);

  profiler_type::stop_pthreading(data->second);

  return NULL;
}

template <class qmci_integrator_type>
void* posix_qmci_integrator<qmci_integrator_type>::start_accumulator_static(void* arg) {
  std::pair<this_type*, int>* data = reinterpret_cast<std::pair<this_type*, int>*>(arg);

  profiler_type::start_pthreading(data->second);

  data->first->start_accumulator(data->second);

  profiler_type::stop_pthreading(data->second);

  return NULL;
}

template <class qmci_integrator_type>
void posix_qmci_integrator<qmci_integrator_type>::start_walker(int id) {
  if (id == 0)
    concurrency << "\n\t\t QMCI starts\n\n";

  walker_type walker(parameters, MOMS, rng_vector[id], id);

  walker.initialize();

  {
    profiler_type profiler("thermalization", "posix-MC-walker", __LINE__, id);
    warm_up(walker, id);
  }

  posix_accumulator_type* acc_ptr(NULL);

  while (acc_finished < nr_accumulators) {
    {
      profiler_type profiler("posix-MC-walker updating", "posix-MC-walker", __LINE__, id);
      walker.do_sweep();
    }

    {
      profiler_type profiler("posix-MC-walker waiting", "posix-MC-walker", __LINE__, id);

      while (acc_finished < nr_accumulators) {
        acc_ptr = NULL;

        {  // checking for available accumulators
          pthread_mutex_lock(&mutex_queue);

          if (!accumulators_queue.empty()) {
            acc_ptr = accumulators_queue.front();
            accumulators_queue.pop();
          }

          pthread_mutex_unlock(&mutex_queue);
        }

        if (acc_ptr != NULL) {
          acc_ptr->update_from(walker);
          acc_ptr = NULL;
          break;
        }

        for (int i = 0; i < parameters.get_additional_steps(); ++i) {
          // std::cout << "\twalker " << id << " is doing some additional steps !!\n";
          profiler_type profiler("additional steps", "posix-MC-walker", __LINE__, id);
          walker.do_step();
        }
      }
    }
  }

  if (id == 0)
    concurrency << "\n\t\t QMCI ends\n\n";
}

template <class qmci_integrator_type>
void posix_qmci_integrator<qmci_integrator_type>::warm_up(walker_type& walker, int id) {
  if (id == 0)
    concurrency << "\n\t\t warm-up starts\n\n";

  for (int i = 0; i < parameters.get_warm_up_sweeps(); i++) {
    walker.do_sweep();

    if (id == 0)
      this->update_shell(i, parameters.get_warm_up_sweeps(), walker.get_configuration().size());
  }

  walker.is_thermalized() = true;

  if (id == 0)
    concurrency << "\n\t\t warm-up ends\n\n";
}

template <class qmci_integrator_type>
void posix_qmci_integrator<qmci_integrator_type>::start_accumulator(int id) {
  posix_accumulator_type accumulator_obj(parameters, MOMS, id);

  accumulator_obj.initialize(DCA_iteration);

  for (int i = 0; i < parameters.get_number_of_measurements(); ++i) {
    pthread_mutex_lock(&mutex_queue);
    accumulators_queue.push(&accumulator_obj);
    pthread_mutex_unlock(&mutex_queue);

    {
      profiler_type profiler("posix-accumulator waiting", "posix-MC-accumulator", __LINE__, id);
      accumulator_obj.wait_for_qmci_walker();
    }

    {
      profiler_type profiler("posix-accumulator accumulating", "posix-MC-accumulator", __LINE__, id);
      if (id == 1)
        this->update_shell(i, parameters.get_number_of_measurements(),
                           accumulator_obj.get_configuration().size());

      accumulator_obj.measure(mutex_queue, accumulators_queue);
    }
  }

  {
    pthread_mutex_lock(&mutex_merge);

    acc_finished++;
    accumulator_obj.sum_to(accumulator);

    pthread_mutex_unlock(&mutex_merge);
  }
}
}

#endif
