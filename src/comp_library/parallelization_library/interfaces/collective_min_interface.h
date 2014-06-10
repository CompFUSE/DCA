//-*-C++-*-

#ifndef COLLECTIVE_MIN_INTERFACE_H
#define COLLECTIVE_MIN_INTERFACE_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  class collective_min_interface
  {
  public:

    collective_min_interface(processor_grouping<LIBRARY>& grouping_ref);
    ~collective_min_interface();

    template<typename scalar_type>
    void min(scalar_type& value);

  private:

    processor_grouping<LIBRARY>& grouping;
  };

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  collective_min_interface<LIBRARY>::collective_min_interface(processor_grouping<LIBRARY>& grouping_ref):
    grouping(grouping_ref)
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  collective_min_interface<LIBRARY>::~collective_min_interface()
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename scalar_type>
  void collective_min_interface<LIBRARY>::min(scalar_type& value)
  {}

}

#endif
