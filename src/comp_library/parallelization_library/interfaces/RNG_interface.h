//-*-C++-*-
//NOTE this file is to be modified to avoid call to std::rand


#ifndef RNG_INTERFACE_H
#define RNG_INTERFACE_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  class RNG_interface
  {
  public:

    RNG_interface(processor_grouping<LIBRARY>& grouping_ref);
    ~RNG_interface();

    void initialize();

    long long get_seed();

    double rand();

  private:

    processor_grouping<LIBRARY>& grouping;

    long long seed;

    double normalization;
  };

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  RNG_interface<LIBRARY>::RNG_interface(processor_grouping<LIBRARY>& grouping_ref):
    grouping(grouping_ref),

    seed(),

    normalization(1)
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  RNG_interface<LIBRARY>::~RNG_interface()
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  void RNG_interface<LIBRARY>::initialize()
  {
    std::srand(3234632461);

    for(int l=0; l<(grouping.get_id()+1)*100; l++)
      seed = std::rand();

    std::srand(seed);

    normalization = 1./double(RAND_MAX);
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  long long RNG_interface<LIBRARY>::get_seed()
  {
    return seed;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  double RNG_interface<LIBRARY>::rand()
  {
    return (std::rand()*normalization);
  }

}

#endif
