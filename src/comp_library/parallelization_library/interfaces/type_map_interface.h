//-*-C++-*-

#ifndef TYPE_MAP_INTERFACE_H
#define TYPE_MAP_INTERFACE_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY, typename scalar_type>
  class type_map_interface
  {
  public:

    static size_t factor() {
      return 1;
    }

    static size_t value() {
      return sizeof(scalar_type);
    }

  };

  /*!
   *  \author Peter Staar
   */
  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY, typename scalar_type>
  class type_map_interface<LIBRARY, std::complex<scalar_type> >
  {
  public:

    static size_t factor() {
      return 2;
    }

    static size_t value() {
      return sizeof(scalar_type);
    }

  };

}
#endif

