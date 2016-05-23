//-*-C++-*-

#ifndef COLLECTIVE_MAX_INTERFACE_MPI_H
#define COLLECTIVE_MAX_INTERFACE_MPI_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<>
  class collective_max_interface<MPI_LIBRARY>
  {
  public:

    collective_max_interface(processor_grouping<MPI_LIBRARY>& grouping_ref);
    ~collective_max_interface();

    template<typename scalar_type>
    void max(scalar_type& value);

  private:

    processor_grouping<MPI_LIBRARY>& grouping;
  };

  collective_max_interface<MPI_LIBRARY>::collective_max_interface(processor_grouping<MPI_LIBRARY>& grouping_ref):
    grouping(grouping_ref)
  {}

  collective_max_interface<MPI_LIBRARY>::~collective_max_interface()
  {}

  template<typename scalar_type>
  void collective_max_interface<MPI_LIBRARY>::max(scalar_type& value)
  {
    scalar_type result;

    MPI_Allreduce(&value,
                  &result,
                  type_map_interface<MPI_LIBRARY, scalar_type>::factor(),
                  type_map_interface<MPI_LIBRARY, scalar_type>::value(),
                  MPI_MAX,
                  grouping.get());

    value = result;
  }

}

#endif
