//-*-C++-*-

#ifndef PROCESSOR_GROUPING_INTERFACE_POSIX_LBRARY_H
#define PROCESSOR_GROUPING_INTERFACE_POSIX_LBRARY_H

namespace COMP_LIB
{
  struct posix_data
  {
  public:

    posix_data();
    ~posix_data();

  public:

    int  id;
    int  nr_threads;

    long seed;

    void* args;
  };

  posix_data::posix_data():
    id(-1),
    nr_threads(-1),
   
    seed(-1),

    args(NULL)
  {}

  posix_data::~posix_data()
  {}

  /*!
   *  \author Peter Staar
   */
  template<>
  class processor_grouping<POSIX_LIBRARY>
  {
  public:

    processor_grouping();
    ~processor_grouping();

    void fork(int N, long SEED, void * (*start_routine)(void *), void *arg);
    void join();

  private:

    std::vector<pthread_t > pthread_vector;
    std::vector<posix_data> data_vector;
  };

  processor_grouping<POSIX_LIBRARY>::processor_grouping():
    pthread_vector(0),
    data_vector   (0)
  {}

  processor_grouping<POSIX_LIBRARY>::~processor_grouping()
  {}

  void processor_grouping<POSIX_LIBRARY>::fork(int N, long seed, void* (*routine)(void *), void *arg)
  {
    srand(seed);

    pthread_vector.resize(N);
    data_vector   .resize(N);

    for(int l=0; l<pthread_vector.size(); l++)
      {
	data_vector[l].id         = l;
	data_vector[l].nr_threads = N;
	
	data_vector[l].seed = rand();

	data_vector[l].args = arg;
	
	pthread_create(&pthread_vector[l], NULL, routine, (void*) (&data_vector[l]));
      }
  }
  
  void processor_grouping<POSIX_LIBRARY>::join()
  {
    for(int l=0; l<pthread_vector.size(); l++)
      pthread_join(pthread_vector[l], NULL);

    pthread_vector.resize(0);
  }

}

#endif
