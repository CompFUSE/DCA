//-*-C++-*-

#ifndef PROCESSOR_GROUPING_INTERFACE_H
#define PROCESSOR_GROUPING_INTERFACE_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  class processor_grouping
  {
  public:

    processor_grouping();
    ~processor_grouping();

    int get_id() ;

    int get_Nr_threads() ;

    int first() ;
    int last () ;

  private:

    int id;         // This processors id within this grouping
    int Nr_threads; // The number of processors within this grouping
  };

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  processor_grouping<LIBRARY>::processor_grouping()
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  processor_grouping<LIBRARY>::~processor_grouping()
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  int processor_grouping<LIBRARY>::get_id()
  {
    return 0;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  int processor_grouping<LIBRARY>::get_Nr_threads()
  {
    return 1;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  int processor_grouping<LIBRARY>::first()
  {
    return 0;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  int processor_grouping<LIBRARY>::last()
  {
    return 0;
  }

}

#endif
