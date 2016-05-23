//-*-C++-*-

#ifndef COLLECTIVE_MIN_INTERFACE_MPI_H
#define COLLECTIVE_MIN_INTERFACE_MPI_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<>
  class collective_min_interface<MPI_LIBRARY>
  {
  public:

    collective_min_interface(processor_grouping<MPI_LIBRARY>& grouping_ref);
    ~collective_min_interface();

    template<typename scalar_type>
    void min(scalar_type& value);

  private:

    processor_grouping<MPI_LIBRARY>& grouping;
  };

  collective_min_interface<MPI_LIBRARY>::collective_min_interface(processor_grouping<MPI_LIBRARY>& grouping_ref):
    grouping(grouping_ref)
  {}

  collective_min_interface<MPI_LIBRARY>::~collective_min_interface()
  {}

  template<typename scalar_type>
  void collective_min_interface<MPI_LIBRARY>::min(scalar_type& value)
  {
    scalar_type result;

    MPI_Allreduce(&value,
                  &result,
                  type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
                  type_map_interface<MPI_LIBRARY, scalar_type>::value(),
                  MPI_MIN,
                  grouping.get());

    value = result;
  }

}

#endif
