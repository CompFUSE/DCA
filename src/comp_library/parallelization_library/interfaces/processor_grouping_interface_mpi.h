//-*-C++-*-

#ifndef PROCESSOR_GROUPING_MPI_H
#define PROCESSOR_GROUPING_MPI_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<>
  class processor_grouping<MPI_LIBRARY>
  {
  public:

    processor_grouping();
    ~processor_grouping();

    void     set();
    MPI_Comm get();

    int get_id();

    int get_Nr_threads();

    int first();
    int last ();

  private:

    int      id;
    int      Nr_threads;

    MPI_Comm MPI_communication;
  };

  processor_grouping<MPI_LIBRARY>::processor_grouping():
    id(-1),
    Nr_threads(-1)
  {}

  processor_grouping<MPI_LIBRARY>::~processor_grouping()
  {}

  void processor_grouping<MPI_LIBRARY>::set()
  {
    MPI_communication = MPI_COMM_WORLD;

    MPI_Comm_size(MPI_COMM_WORLD, &Nr_threads);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
  }

  MPI_Comm processor_grouping<MPI_LIBRARY>::get()
  {
    return MPI_communication;
  }

  int processor_grouping<MPI_LIBRARY>::get_id()
  {
    assert(id>-1);
    return id;
  }

  int processor_grouping<MPI_LIBRARY>::get_Nr_threads()
  {
    assert(Nr_threads>-1);
    return Nr_threads;
  }

  int processor_grouping<MPI_LIBRARY>::first()
  {
    return 0;
  }

  int processor_grouping<MPI_LIBRARY>::last ()
  {
    return Nr_threads-1;
  }

}

#endif
