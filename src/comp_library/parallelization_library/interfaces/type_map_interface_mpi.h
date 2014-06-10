//-*-C++-*-
#ifndef TYPE_MAP_INTERFACE_MPI_H
#define TYPE_MAP_INTERFACE_MPI_H

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<typename scalar_type>
  class type_map_interface<MPI_LIBRARY, scalar_type>
  {
  public:

    static size_t factor() {
      return 1;
    }

    static MPI_Datatype value() {
      return  MPI_PACKED;
    }
  };

  /*!
   *  \author Peter Staar
   */
  template<>
  class type_map_interface<MPI_LIBRARY, bool>
  {
  public:

    static size_t factor() {
      return 1;
    }

    static MPI_Datatype value() {
      return  MPI_INT;
    }
  };

  /*!
   *  \author Peter Staar
   */
  template<>
  class type_map_interface<MPI_LIBRARY, char>
  {
  public:

    static size_t factor() {
      return 1;
    }

    static MPI_Datatype value() {
      return  MPI_CHAR;
    }
  };

  /*!
   *  \author Peter Staar
   */
  template<>
  class type_map_interface<MPI_LIBRARY, int>
  {
  public:

    static size_t factor() {
      return 1;
    }

    static MPI_Datatype value() {
      return  MPI_INT;
    }
  };

  /*!
   *  \author Peter Staar
   */
  template<>
  class type_map_interface<MPI_LIBRARY, size_t>
  {
  public:

    static size_t factor() {
      return 1;
    }

    static MPI_Datatype value() {
      return  MPI_UNSIGNED_LONG;
    }
  };

  /*!
   *  \author Peter Staar
   */
  template<>
  class type_map_interface<MPI_LIBRARY, float>
  {
  public:

    static size_t factor() {
      return 1;
    }

    static MPI_Datatype value() {
      return  MPI_FLOAT;
    }
  };

  /*!
   *  \author Peter Staar
   */
  template<>
  class type_map_interface<MPI_LIBRARY, double>
  {
  public:

    static size_t factor() {
      return 1;
    }

    static MPI_Datatype value() {
      return  MPI_DOUBLE;
    }
  };

  /*!
   *  \author Peter Staar
   */
  template<>
  class type_map_interface<MPI_LIBRARY, std::complex<float> >
  {
  public:

    static size_t factor() {
      return 2;
    }

    static MPI_Datatype value() {
      return  MPI_FLOAT;
    }
  };

  /*!
   *  \author Peter Staar
   */
  template<>
  class type_map_interface<MPI_LIBRARY, std::complex<double> >
  {
  public:

    static size_t factor() {
      return 2;
    }

    static MPI_Datatype value() {
      return  MPI_DOUBLE;
    }
  };

}

#endif
