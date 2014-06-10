//-*-C++-*-

#ifndef PACKING_INTERFACE_HEADER
#define PACKING_INTERFACE_HEADER

namespace COMP_LIB
{
  /*!
   *  \author Peter Staar
   */
  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  class packing_interface
  {
  public:

    packing_interface(processor_grouping<LIBRARY>& grouping_ref);
    ~packing_interface();

    /************************************
     ***  size
     ************************************/

    template<typename T>
    size_t get_buffer_size(T& item);

    size_t get_buffer_size(std::string str);

    template<typename T>
    size_t get_buffer_size(std::vector<T>& v);

    template<typename T, class dmn_type>
    size_t get_buffer_size(function<T, dmn_type>& f);

    /************************************
     ***  pack_unpack
     ************************************/

    template<typename T>
    void pack_unpack(bool packing, int* buffer, int size, int& off_set, T& item);

  private:

    /************************************
     ***  pack
     ************************************/

    template<typename T>
    void pack(int* buffer, int size, int& off_set, T& item );

    void pack(int* buffer, int size, int& off_set, std::string& str );

    template<typename T>
    void pack(int* buffer, int size, int& off_set, std::vector<T>& v);

    template<typename T, class dmn_type>
    void pack(int* buffer, int size, int& off_set, function<T, dmn_type>& f);

    /************************************
     ***  unpack
     ************************************/

    template<typename T>
    void unpack(int* buffer, int size, int& off_set, T& item);

    void unpack(int* buffer, int size, int& off_set, std::string& str);

    template<typename T>
    void unpack(int* buffer, int size, int& off_set, std::vector<T>& v);

    template<typename T, class dmn_type>
    void unpack(int* buffer, int size, int& off_set, function<T, dmn_type>& f);

  private:

    processor_grouping<LIBRARY>& grouping;
  };

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  packing_interface<LIBRARY>::packing_interface(processor_grouping<LIBRARY>& grouping_ref):
    grouping(grouping_ref)
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  packing_interface<LIBRARY>::~packing_interface()
  {}

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T>
  size_t packing_interface<LIBRARY>::get_buffer_size(T& item)
  {
    return sizeof(T);
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  size_t packing_interface<LIBRARY>::get_buffer_size(std::string str)
  {
    return str.size()*sizeof(char);
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T>
  size_t packing_interface<LIBRARY>::get_buffer_size(std::vector<T>& v)
  {
    return v.size()*size(v[0]);
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T, class dmn_type>
  size_t packing_interface<LIBRARY>::get_buffer_size(function<T, dmn_type>& f)
  {
    return f.size()*size(f(0));
  }

  /************************************
   ***  pack_unpack
   ************************************/

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T>
  void packing_interface<LIBRARY>::pack_unpack(bool packing, int* buffer, int size, int& off_set, T& item)
  {
    if(packing)
      pack  (buffer, size, off_set, item);
    else
      unpack(buffer, size, off_set, item);
  }

  /************************************
   ***  pack
   ************************************/

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T>
  void packing_interface<LIBRARY>::pack(int* buffer, int size, int& off_set, T& item)
  {
    //   size_t object_size = size(item);
    //   int* ptr = static_cast<int*>(&(item));
    //   memcpy(buffer+offset, &ptr, object_size);
    //   off_set += object_size;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  void packing_interface<LIBRARY>::pack(int* buffer, int size, int& off_set, std::string& str)
  {
    //   for(size_t l=0; l<str.size(); l++)
    //     pack(buffer, size, off_set, str[l])
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T>
  void packing_interface<LIBRARY>::pack(int* buffer, int size, int& off_set, std::vector<T>& v)
  {
    //   for(size_t l=0; l<v.size(); l++)
    //     pack(buffer, size, off_set, v[l])
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T, class dmn_type>
  void packing_interface<LIBRARY>::pack(int* buffer, int size, int& off_set, function<T, dmn_type>& f)
  {
    //   for(int l=0; l<f.size(); l++)
    //     pack(buffer, size, off_set, f(l));
  }

  /************************************
   ***  unpack
   ************************************/

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T>
  void packing_interface<LIBRARY>::unpack(int* buffer, int size, int& off_set, T& item)
  {
    //   size_t object_size = size(item);
    //   int* ptr = static_cast<int*>(&(item));
    //   memcpy(&ptr, buffer+offset, object_size);
    //   off_set += object_size;
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  void packing_interface<LIBRARY>::unpack(int* buffer, int size, int& off_set, std::string& str)
  {
    //   for(size_t l=0; l<str.size(); l++)
    //     unpack(buffer, size, off_set, str[l])
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T>
  void packing_interface<LIBRARY>::unpack(int* buffer, int size, int& off_set, std::vector<T>& v)
  {
    //   for(size_t l=0; l<v.size(); l++)
    //     unpack(buffer, size, off_set, v[l])
  }

  template<PARALLELIZATION_LIBRARY_NAMES LIBRARY>
  template<typename T, class dmn_type>
  void packing_interface<LIBRARY>::unpack(int* buffer, int size, int& off_set, function<T, dmn_type>& f)
  {
    //   for(int l=0; l<f.size(); l++)
    //     unpack(buffer, size, off_set, f(l));
  }

}

#endif
