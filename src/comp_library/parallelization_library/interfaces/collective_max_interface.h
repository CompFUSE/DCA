//-*-C++-*-

#ifndef COLLECTIVE_MAX_INTERFACE_H
#define COLLECTIVE_MAX_INTERFACE_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  class collective_max_interface
  {
  public:

    collective_max_interface(processor_grouping<LIBRARY>& grouping_ref);
    ~collective_max_interface();

    template<typename scalar_type>
    void max(scalar_type& value);

  private:

    processor_grouping<LIBRARY>& grouping;
  };

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  collective_max_interface<LIBRARY>::collective_max_interface(processor_grouping<LIBRARY>& grouping_ref):
    grouping(grouping_ref)
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  collective_max_interface<LIBRARY>::~collective_max_interface()
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename scalar_type>
  void collective_max_interface<LIBRARY>::max(scalar_type& value)
  {}

}

#endif
