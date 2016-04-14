//-*-C++-*-

#ifndef PARALLELIZATION_LIBRARY_TEMPLATE_H
#define PARALLELIZATION_LIBRARY_TEMPLATE_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  class parallelization : public print_on_shell_interface<LIBRARY>,
                          public packing_interface       <LIBRARY>,
                          public collective_sum_interface<LIBRARY>
  {

  public:

    parallelization(int argc, char *argv[]);
    ~parallelization();

    int id();

    int number_of_processors();

    int first();
    int last ();

    template<typename object_type>
    bool broadcast(object_type& object, int root_id=0);

    template<typename object_type>
    bool broadcast_object(object_type& object, int root_id=0);

    template<typename domain_t>
    std::pair<int, int> get_bounds(domain_t& dmn);

  private:

    processor_grouping<LIBRARY> group;
  };

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  parallelization<LIBRARY>::parallelization(int /*argc*/, char */*argv*/[]):
    print_on_shell_interface<LIBRARY>(group),
    packing_interface       <LIBRARY>(group),
    collective_sum_interface<LIBRARY>(group)
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  parallelization<LIBRARY>::~parallelization()
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  int parallelization<LIBRARY>::id()
  {
    return 0;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  int parallelization<LIBRARY>::number_of_processors()
  {
    return 1;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  int parallelization<LIBRARY>::first()
  {
    return 0;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  int parallelization<LIBRARY>::last()
  {
    return 0;
  }


  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename object_type>
  bool parallelization<LIBRARY>::broadcast(object_type& /*object*/, int /*root_id*/)
  {
    return true;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename object_type>
  bool parallelization<LIBRARY>::broadcast_object(object_type& /*object*/, int /*root_id*/)
  {
    return true;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename domain_t>
  std::pair<int,int> parallelization<LIBRARY>::get_bounds(domain_t& dmn)
  {
    int size = dmn.get_size();

    std::pair<int, int> bounds;
    bounds.first  = 0;
    bounds.second = size;

    return bounds;
  }

}

#endif
