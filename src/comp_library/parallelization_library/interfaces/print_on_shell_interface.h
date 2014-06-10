//-*-C++-*-

#ifndef PRINT_ON_SHELL_INTERFACE_H
#define PRINT_ON_SHELL_INTERFACE_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  class print_on_shell_interface
  {
  public:

    print_on_shell_interface(processor_grouping<LIBRARY>& grouping_ref);
    ~print_on_shell_interface();

    template<class whatever>
    void operator<<(whatever& whtvr);

    template<class whatever>
    void operator<<(const whatever& whtvr);

    template<class whatever>
    void operator<<(whatever* whtvr);

    template<class whatever>
    void operator<<(const whatever* whtvr);

  private:

    processor_grouping<LIBRARY>& grouping;
  };

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  print_on_shell_interface<LIBRARY>::print_on_shell_interface(processor_grouping<LIBRARY>& grouping_ref):
    grouping(grouping_ref)
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  print_on_shell_interface<LIBRARY>::~print_on_shell_interface()
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<class whatever>
  void print_on_shell_interface<LIBRARY>::operator<<(whatever& whtvr)
  {
    if(grouping.get_id() == grouping.first())
      cout << whtvr;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<class whatever>
  void print_on_shell_interface<LIBRARY>::operator<<(const whatever& whtvr)
  {
    if(grouping.get_id() == grouping.first())
      cout << whtvr;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<class whatever>
  void print_on_shell_interface<LIBRARY>::operator<<(whatever* whtvr)
  {
    if(grouping.get_id() == grouping.first())
      cout << whtvr;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<class whatever>
  void print_on_shell_interface<LIBRARY>::operator<<(const whatever* whtvr)
  {
    if(grouping.get_id() == grouping.first())
      cout << whtvr;
  }

}

#endif
